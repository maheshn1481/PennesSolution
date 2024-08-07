
tic
clear all
clc
close all
% Parameters
Rt = 0.005;                 % Radius of tumor in which MNP heating is present
R = 5*Rt;                   % Radius of sphere (m)
dr = 0.0005;                 % Spatial element size
r = 0:dr:R;                 % Spatial locations
N = length(r);              % Number of nodes
dt = 1;                     % Time step (s)
t_final = 600;              % Total simulation time (s)
time = 0:dt:t_final;        % Simulation time points
t_points = length(time);     % Number of time steps 

% Tissue properties
rho = 1000;        % Density (kg/m^3)
cp = 3700;          % Specific heat (J/kg·K)
k = 0.5;           % Thermal conductivity (W/m·K)
w_b = 0.0005;      % Blood perfusion rate (s^-1)
rho_b = 1000;      % Blood density (kg/m^3)
c_b = 3600;        % Blood specific heat (J/kg·K)
T_a = 37;          % Arterial blood temperature (°C)
%q_m = 420000;         % Metabolic heat generation (W/m^3)
mnp_loc = r<=Rt;
q_m = 500000*mnp_loc;      % Metabolic heat generation (W/m^3)

% Boundary condition
h = 10;            % Convective heat transfer coefficient (W/m^2·K)
T_inf = 25;        % Ambient temperature (°C)


% Initialize solution array
T_all = zeros(N, t_points);

% Initial condition
T = ones(N, 1) * 37;  % Initial temperature (°C)
T_all(:, 1) = T;

% Time stepping
for step = 2:t_points
    % Set up tridiagonal system
    a = zeros(N, 1);
    b = zeros(N, 1);
    c = zeros(N, 1);
    d = zeros(N, 1);

    % Interior points
    for i = 2:N-1
        a(i) = -(k/rho*cp)*(dt/dr^2)+(k/rho*cp)*(dt/(r(i)*dr));
        b(i) = 1 + 2*(k/rho*cp)*(dt/dr^2) + (w_b*rho_b*c_b*dt)/(rho*cp);
        c(i) = -(k/rho*cp)*(dt/dr^2)-(k/rho*cp)*(dt/(r(i)*dr));
        % d(i) = T(i) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m*dt)/(rho*cp)+(q_mnp(i)*dt)/(rho*cp);
        d(i) = T(i) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m(i)*dt)/(rho*cp);
    end

    % Boundary conditions
    % r = 0 (symmetry)
    b(1) = 1 + 2*(k/rho*cp)*(dt/dr^2) + (w_b*rho_b*c_b*dt)/(rho*cp);
    c(1) = -2*(k/rho*cp)*(dt/dr^2);
    d(1) = T(1) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m*dt)/(rho*cp)+(q_mnp(1)*dt)/(rho*cp);
    % d(1) = T(1) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m(1)*dt)/(rho*cp);

    % % r = R (convective)
    % a(N) = -k*dt/(rho*cp*dr^2);
    % b(N) = 1 + k*dt/(rho*cp*dr^2) + h*dt/(rho*cp*dr) + w_b*rho_b*c_b*dt/(rho*cp);
    % d(N) = T(N) + w_b*rho_b*c_b*dt/(rho*cp)*T_a + q_m*dt/(rho*cp) + h*dt/(rho*cp*dr)*T_inf+(q_mnp(N)*dt)/(rho*cp);
    
    % r = R (No flux)
    % a(N) = 1 + 2*(k/rho*cp)*(dt/dr^2) + (w_b*rho_b*c_b*dt)/(rho*cp);
    % b(N) = -2*(k/rho*cp)*(dt/dr^2);
    % % d(N) = T(N) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m*dt)/(rho*cp)+(q_mnp(N)*dt)/(rho*cp);
    % d(N) = T(N) + (w_b*rho_b*c_b*dt*T_a)/(rho*cp) + (q_m(N)*dt)/(rho*cp);

    % % r = R (Dirichlet)
    % a(N) = 1;
    % d(N) = T_a;

    % Solve tridiagonal system
    T = tridiag(a, b, c, d);

    % Store solution
    T_all(:, step) = T;
end

% Plot results
figure;

% Final temperature distribution
subplot(2,1,1);
plot(r, T_all(:, end));
xlabel('Radius (m)');
ylabel('Temperature (°C)');
title('Final Temperature Distribution in Sphere');
grid on;

% Temperature evolution at different radii
subplot(2,1,2);
plot(time, T_all(1,:), time, T_all(round(N/2),:), time, T_all(end,:));
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Temperature Evolution');
legend('Center', 'Midpoint', 'Surface');
grid on;

% Create animation
% figure;
% for step = 1:t_points
%     plot(r, T_all(:, step));
%     xlabel('Radius (m)');
%     ylabel('Temperature (°C)');
%     title(sprintf('Temperature Distribution at t = %.2f s', time(step)));
%     ylim([min(T_all(:)) max(T_all(:))]);
%     grid on;
%     drawnow;
%     pause(0.01);
% end
% end
toc

function x = tridiag(a, b, c, d)
% Solves a tridiagonal system Ax = d where A has diagonals a, b, c
n = length(d);
x = zeros(n, 1);

% Forward elimination
for i = 2:n
    m = a(i) / b(i-1);
    b(i) = b(i) - m * c(i-1);
    d(i) = d(i) - m * d(i-1);
end

% Back substitution
x(n) = d(n) / b(n);
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end
end

