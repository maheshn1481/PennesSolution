tic
clc
clear
close all
% MATLAB code to solve the one-dimensional transient Pennes bioheat equation
% in spherical coordinates using an implicit time discretization scheme.

% Parameters
rho1 = 1660; % density of region 1 (kg/m^3)
c1 = 2540; % specific heat of region 1 (J/kg*K)
k1 = 0.778; % thermal conductivity of region 1 (W/m*K)
wb1 = 0.009; % blood perfusion rate of region 1 (1/s)
qm1 = 29000; % metabolic heat generation (W/m^3)
qe1 = 6.15e6; % external heat source (W/m^3)

rho2 = 1000; % density of region 2 (kg/m^3)
c2 = 3720; % specific heat of region 2 (J/kg*K)
k2 = 0.642; % thermal conductivity of region 2 (W/m*K)
wb2 = 0.0018; % blood perfusion rate of region 2 (1/s)
qm2 = 450; % metabolic heat generation (W/m^3)
qe2 = 0; % external heat source (W/m^3)

rhoi = 2*rho1*rho2/(rho1+rho2);
ci = 2*c1*c2/(c1+c2);
ki = 2*k1*k2/(k1+k2);
wbi = 2*wb1*wb2/(wb1+wb2);
qmi = 2*qm1*qm2/(qm1+qm2);


rhob = 1000; % blood density (kg/m^3)
cb = 4180; % specific heat of blood (J/kg*K)
Tb = 37; % arterial blood temperature (C)



r_i = 0.00315; % interface radius (m)
R = 0.015; % outer radius (m)
T_s = 37; % specified surface temperature (C)
T0 = 37; % initial temperature (C)

dr = r_i/100; % spatial step size (m)
dt = 0.1; % time step size (s)
t_end = 1000; % end time (s)

% Discretization
r = 0:dr:R;
N = length(r);
M = t_end / dt;



% Find interface index
i_interface = find(r_i==r);

alpha1 = (k1*dt)/(rho1*c1*dr^2);
alpha2 = (k2*dt)/(rho2*c2*dr^2);
beta1 = rhob*cb*wb1*dt/(rho1*c1);
beta2 = rhob*cb*wb2*dt/(rho2*c2);
alphai = (ki*dt)/(rhoi*ci*dr^2);
betai = rhob*cb*wbi*dt/(rhoi*ci);
Lower = zeros(N, 1);
Main  = zeros(N, 1);
Upper = zeros(N, 1);
for i = 1:N
    if i==1
        Lower(i) = 0;
        Main(i) = 1 + 2 * alpha1 + beta1;
        Upper(i) = -2*alpha1;
    elseif i~=1 && i<i_interface
        Lower(i) = -alpha1 * (1 - dr / r(i));
        Main(i) = 1 + 2 * alpha1 + beta1;
        Upper(i) = -alpha1 * (1 + dr / r(i));
    elseif i==i_interface
        Lower(i) = -k1;
        Main(i) = k1+k2;
        Upper(i) = -k2;
    elseif i~=N && i>i_interface
        Lower(i) = -alpha2 * (1 - dr / r(i));
        Main(i) = 1 + 2 * alpha2 + beta2;
        Upper(i) = -alpha2 * (1 + dr / r(i));
    elseif i==N
        Lower(i) = 0;
        Main(i) = 1;
        Upper(i) = 0;
    end
end
% Initialize temperature arrays
T = T0 * ones(N, M+1);  % Store all time steps, including initial condition

for t = 1:M
    Force = zeros(N,1);
    for i = 1:N
        if i==1
            Force(i) = T(i, t) + dt * (beta1 * Tb + (qm1 + qe1) / (rho1 * c1));
        elseif i~=1 && i<i_interface
            Force(i) = T(i, t) + dt * (beta1 * Tb + (qm1 + qe1) / (rho1 * c1));
        elseif i==i_interface
            Force(i) = 0;
        elseif i~=N && i>i_interface
            Force(i) = T(i, t) + dt * (beta2 * Tb + (qm2 + qe2) / (rho2 * c2));
        elseif i==N
            Force(i) = T_s;
        end
    end
    T(:,t+1) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
end

T = T-37;

% Plot the temperature distribution at different times
figure;
times_to_plot = [1, round(M/4), round(M/2), round(3*M/4), M+1];
colors = {'b', 'g', 'r', 'c', 'm'};
legends = cell(1, length(times_to_plot));

for i = 1:length(times_to_plot)
    t_index = times_to_plot(i);
    plot(r, T(:, t_index), 'Color', colors{i}, 'LineWidth', 2);
    hold on;
    legends{i} = sprintf('t = %.2f s', (t_index-1)*dt);
end

xlabel('Radius (m)');
ylabel('Temperature (C)');
title('Temperature Distribution in Spherical Coordinates');
grid on;

% Add a vertical line to show the interface
line([r_i r_i], ylim, 'Color', 'k', 'LineStyle', '--');
legend([legends, 'Interface'], 'Location', 'best');

% Create a surface plot to show temperature evolution over time and space
figure;
[R, Time] = meshgrid(r, 0:dt:t_end);
surf(R, Time, T');
xlabel('Radius (m)');
ylabel('Time (s)');
zlabel('Temperature (C)');
title('Temperature Evolution in Space and Time');
colorbar;
shading interp;

% Add a plane to show the interface
hold on;
[X, Y] = meshgrid(linspace(0, t_end, 100), linspace(min(T(:)), max(T(:)), 100));
Z = r_i * ones(size(X));
surf(Z, X, Y, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold off;
toc

%% Tri-Diagonal Matrix Alogorithm for Linear System of Equation Solver
function x = thomas_algorithm(Lower, Main, Upper, Force)
n = length(Force);
c_star = zeros(n, 1);
d_star = zeros(n, 1);
x = zeros(n, 1);

c_star(1) = Upper(1) / Main(1);
d_star(1) = Force(1) / Main(1);

for i = 2:n
    c_star(i) = Upper(i) / (Main(i) - Lower(i) * c_star(i-1));
    d_star(i) = (Force(i) - Lower(i) * d_star(i-1)) / (Main(i) - Lower(i) * c_star(i-1));
end

x(n) = d_star(n);
for i = n-1:-1:1
    x(i) = d_star(i) - c_star(i) * x(i+1);
end
end