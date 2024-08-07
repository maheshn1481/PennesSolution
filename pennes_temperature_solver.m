function [T,r,t,q_mnp] = pennes_temperature_solver(modelInput,heat_source,Blood,Tumor,Mag_Field)
%%% Domain Parameters
RT = modelInput.tumor_size;% [m] Radius of the Tumor
R = 5*RT;% [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;% [m3] Tumor Volume
dr = RT/50;% [m] Spatial discretization step size
r = (0:dr:R)';% [m] Spatial loactions
N = length(r);% [-] Number of spatial nodes
t_loc = r<=RT;% [-] Tumor r indices

%%% Time discretization: Implicit method is employed
dt = 1;% Time step [s]
t_end = modelInput.treatment_time;% End time [s]
t = 0:dt:t_end;% Discrete times
TS = length(t);% Number of time steps

%%% Initial Conditions
T_initial = modelInput.initial_condition;% [°C] Initial Condition Temperature

%%% Blood Properties
rho_b = Blood.rho;% [kg/m3] Blood Density
cp_b = Blood.cp;% [J/kg/K] Blood Specific Heat
T_b = Blood.T;% [°C] Blood Arterial Temperature
T_bK = T_b + 273.15;% [K] Blood Arterial Temperature

%%% Tumor Properties
rho_t = Tumor.rho;               % [kg/m3] Tumor Density
cp_t = Tumor.cp;                % [J/kg/K] Tumor Specific Heat
k_t = Tumor.k;                % [W/m/K] Tumor Thermal Conductivity
w_t = Tumor.w;              % [1/s] Tumor Blood Perfusion Rate
Qm_t = Tumor.Qm;                   % [W/m3] Tumor Metabolic Heat
alpha_t = k_t/(rho_t*cp_t); % Thermal diffusivity [m^2/s]



%%% Magnetic Field Parameters
H = Mag_Field.Amplitude;  % [A/m] Magnetic Field Amplitude H = 10 kA/m
f = Mag_Field.Frequency;    % [Hz] Magnetic Field Frequency f = 100 kHz
Hf = H*f;   % [A.Hz/m] Amplitude and Frequency product

%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 1.6e-8;                 % [m] MNP Diameter d = 16 nm
delta_l = 1e-9;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
rho_mnp = 5180;                 % [kg/m3] MNP Density
cp_mnp = 670;                   % [J/kg/K] MNP Specific Heat
Md = 4.46e5;                    % [A/m] Domain Magnetization Value Md = 446 kA/m
Keff = 20000;                   % [J/m3] Effective Anisotropy Constant Keff = 20 kJ/m3

%%% MNP Dose and Magnetic Fluid Parameters
MNP_dose = 4;                   % [kg/m3] MNP Dose MNP mass used per unit vol of tumor = 4 mg/cm3; % Conversion 1 mg/cm3 = 1 kg/m3
MNP_mass = MNP_dose*vol_t;      % [kg] MNP mass injected into the tumor; Varies with MNP Dose
MNP_vol = MNP_mass/rho_mnp;     % [m3] volume of MNP in the Magnetic Fluid
MF_inj_vol = 1;                 % [ml, milliliters] Magnetic Fluid Volume injected into the Tumor
MF_vol = MF_inj_vol*1e-6;       % [m3] Magnetic Fluid Volume Injected in to the Tumor % Converion Factor for ml to m3 is 1e-6
rho_cf = 1000;                  % [kg/m3] Density of the Carrier Fluid
mu_cf = 1e-3;                   % [Pa.s] Viscocity of the Carrier Fluid
MNP_vol_frac = MNP_vol/MF_vol;  %[-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7;              % [H/m] Permeability of Free Space
kB = 1.38e-23;              % [J/K] Blotzmann Constant
MNP_svol = pi*(d_mnp^3)/6;  % [m3] Solid Volume of Individual Mangetic Nanoparticle
MNP_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9;                % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f;             % [rad/s] Angular frequency of applied magnetic field
alpha_CF = 0.55;            % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Keff*MNP_svol/(kB*T_bK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*MNP_hvol)/(kB*T_bK);                % [s] Brownian Relaxation Time of MNP
tauN = sqrt(pi)*tau0*exp(gamma)/(2*sqrt(gamma));    % [s] Neel Relaxation Time of MNP
tauE = tauB*tauN/(tauB+tauN);                       % [s] Effective Relaxation Time of MNP
zeta = mu0*Md*H*MNP_svol/(kB*T_bK);                 % [-] Langevin parameter
X_i = mu0*(Md^2)*MNP_vol_frac*MNP_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
SAR = P/(rho_mnp*MNP_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
SAR_grams = SAR/1000;                               % [W/g] MNP Specific Absorption Rate

%%% Heat Source by MNP in Tumor
MNP_conc = MNP_mass/vol_t;      % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
Q_MNP = alpha_CF*MNP_conc*SAR;  % [W/m3] Heat Generation by MNP in Tissue

%% Implicit Finite Difference Scheme
T = zeros(N, TS);
T(:,1) = T_initial; % Initial Condition
%%% Coefficient Matrix Three-Diagonal Generation
Lower = zeros(N, 1);
Main  = zeros(N, 1);
Upper = zeros(N, 1);
for i = 1:N
    if i ==1 % Symmetry Node
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rho_t*cp_t);
        Upper(i) = -2*alpha_t*dt/(dr^2);
    elseif i == N % Fixed Temperature Node
        Lower(i) = -0;
        Main(i) = 1;
    else % Internal Nodes
        Lower(i) = -alpha_t*dt/(dr^2) + alpha_t*dt/(r(i)*dr);
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rho_t*cp_t);
        Upper(i) = -alpha_t*dt/(dr^2) - alpha_t*dt/(r(i)*dr);
    end
end
%% Solver
switch heat_source.type
    case 'Constant'
        %%% Assign Q_MNP to the center region only
        q_mnp = Q_MNP*t_loc;            % MNP heat generation [W/m^3] On for all simulation time
        %%% Time Iterations
        for n = 2:TS
            Force = zeros(N,1); % Force vector
            Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1)*dt)/(rho_t*cp_t);
            Force(N) = T_b;
            T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
        end
        % %%% Plot the results
        % figure('Position', [100, 100, 800, 400])
        % subplot(1,2,1)
        % plot_radius = [0,RT,R];
        % legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
        % plot_rad_idx = (plot_radius/dr)+1;
        % plot(t,T(plot_rad_idx,:), 'LineWidth',2)
        % xlabel('Time, t [s]');
        % ylabel('Temperature, T [°C]');
        % xlim([0,t_end])
        % xticks(0:200:t_end);
        % legend(legend_String, Location='south')
        % title('Temperature Elevations');
        % grid on;
        % subplot(1,2,2)
        % plot_time = [0,60,120,t_end];
        % legend_String = string(plot_time)+[" s"," s"," s"," s"];
        % plot_ind = (plot_time/dt)+1;
        % plot(r,T(:,plot_ind),'LineWidth',2)
        % hold on
        % xlabel('Radial distance, r [m]');
        % xlim([0,R])
        % xticks(0:R/5:R);
        % ylabel('Temperature, T [°C]');
        % legend(legend_String)
        % title('Spatial Temperature profile');
    case 'Pulsed'
        %%% Creation of Pulsating Function in time to Turn ON/OFF MNP Heat source
        cycle_on = heat_source.on_time;      % [s] Cycle ON time
        t_vec = t;          % [s] Time vector
        step_height = 1;    % [-] Pulse Amplitude
        y = step_height * mod(floor(t_vec / cycle_on), 2);
        y = step_height - y; % Invert to start from peak (1)
        y(end) =[];
        y = [1,y]; % Pulsating Profile
        heating_map = t_loc*y; % Used to assigne the Q_MNP within tumor and at ON times
        q_mnp = Q_MNP*heating_map;        % MNP heat generation [W/m^3] On/Off Pulsating with time

        %%% Time Iterations
        for n = 2:TS
            Force = zeros(N,1); % Force vector
            Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1,n)*dt)/(rho_t*cp_t);
            Force(N) = T_b;
            T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
        end
        % %%% Plot the results
        % figure('Position', [100, 100, 1200, 400])
        % subplot(1,3,1)
        % plot(t_vec, y*Q_MNP, 'LineWidth', 2);
        % xlabel('Time [s]');
        % ylabel('Qmnp [W/m^3]');
        % title('Pulsed Heat Source');
        % grid on;
        % xticks(0:2*cycle_on:t_end);
        % xlim([0 4*cycle_on]);
        % subplot(1,3,2)
        % plot_radius = [0,RT,R];
        % legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
        % plot_rad_idx = (plot_radius/dr)+1;
        % plot(t,T(plot_rad_idx,:), 'LineWidth',2)
        % xlabel('Time, t [s]');
        % ylabel('Temperature, T [°C]');
        % xlim([0,t_end])
        % xticks(0:200:t_end);
        % legend(legend_String, Location='south')
        % title('Temperature Elevations');
        % grid on;
        % subplot(1,3,3)
        % plot_time = [0,60,120,t_end];
        % legend_String = string(plot_time)+[" s"," s"," s"," s"];
        % plot_ind = (plot_time/dt)+1;
        % plot(r,T(:,plot_ind),'LineWidth',2)
        % hold on
        % xlabel('Radial distance, r [m]');
        % xlim([0,R])
        % xticks(0:R/5:R);
        % ylabel('Temperature, T [°C]');
        % legend(legend_String)
        % title('Spatial Temperature profile');
        % grid on;
end
end
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

