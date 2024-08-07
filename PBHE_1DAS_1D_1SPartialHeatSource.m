tic
clear; close all; clc;
%%% Uses TDMA algorithm to solve Discretized Linear Equations %%%
%% Problem Parameters

%%% Domain Parameters
RT = 0.005;                 % [m] Radius of the Tumor
R = 5*RT; % [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;    % [m3] Tumor Volume
dr = RT/50;                % [m] Spatial discretization step size
r = (0:dr:R)';            % [m] Spatial loactions
N = length(r);            % [-] Number of spatial nodes
t_loc = r<=RT;            % [-] Tumor r indices



%%% Time discretization: Implicit method is employed
dt = 2;             % Time step [s]
t_end = 1200;        % End time [s]
t = 0:dt:t_end;     % Discrete times
TS = length(t);     % Number of time steps


%%% Blood Properties
rho_b = 1000;           % [kg/m3] Blood Density
cp_b = 3840;            % [J/kg/K] Blood Specific Heat
k_b = 0.50;             % [W/m/K] Blood Thermal Conductivity
T_b = 37;               % [°C] Blood Arterial Temperature
T_bK = T_b + 273.15;    % [K] Blood Arterial Temperature


%%% Tumor Properties
rho_t = 1045;   % [kg/m3] Tumor Density
cp_t = 3600;    % [J/kg/K] Tumor Specific Heat
k_t = 0.527;     % [W/m/K] Tumor Thermal Conductivity
w_t = 9/rho_b;  % [1/s] Tumor Blood Perfusion Rate
Qm_t = 0;    % [W/m3] Tumor Metabolic Heat
alpha_t = k_t/(rho_t*cp_t); % Thermal diffusivity [m^2/s]


%%% Magnetic Field Parameters
H = 10000;  % [A/m] Magnetic Field Amplitude H = 10 kA/m
f = 1e5;    % [Hz] Magnetic Field Frequency f = 100 kHz
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
%%% Assign Q_MNP to the tumor grids only
q_mnp = Q_MNP*t_loc;            % MNPc heat generation [W/m^3]

%% Implicit Finite Difference Scheme
T = zeros(N, TS);
T(:,1) = T_b;
Lower = zeros(N, 1);
Main  = zeros(N, 1);
Upper = zeros(N, 1);
Force = zeros(N, 1);
% %%% Considering Symmetry at r = 0 and fixed Temperature at r = R
for i = 1:N
    if i ==1
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rho_t*cp_t);
        Upper(i) = -2*alpha_t*dt/(dr^2);
    elseif i == N
        Lower(i) = -0;
        Main(i) = 1;
    else
        Lower(i) = -alpha_t*dt/(dr^2) + alpha_t*dt/(r(i)*dr);
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rho_t*cp_t);
        Upper(i) = -alpha_t*dt/(dr^2) - alpha_t*dt/(r(i)*dr);
    end
end

for n = 2:TS
    Force = zeros(N,1);
    Force(1:N-1)  = T(1:N-1,n-1) + rho_b*cp_b.*w_t*dt*T_b/(rho_t*cp_t) + Qm_t*dt/(rho_t*cp_t) + q_mnp(1:N-1)*dt/(rho_t*cp_t);
    Force(N) = T_b;
    T(:,n) = tridiag(Lower, Main, Upper, Force);
end

%% Plot the results
close all
figure
subplot(1,2,1)
plot(t,T(1,:),t,T(find(RT==r),:), 'LineWidth',2)
xlabel('Time, t [s]');
ylabel('Temperature, T [°C]');
xlim([0,t_end])
xticks(0:200:t_end);
title('Temperature Elevations');
% legend('0','0.25R','0.5R','0.75R','R',Location='best')
grid on;
subplot(1,2,2)
plot(r,T(:,end),'LineWidth',2)
xlabel('Radial distance, r [m]');
xlim([0,R])
xticks(0:R/5:R);
ylabel('Temperature, T [°C]');
title('Spatial Temperature profile');
% legend('0 s', '10 s','20 s','1200 s')
%%
toc
%%
function x = tridiag(Lower, Main, Upper, Force)
% Solves a tridiagonal system Ax = Force where A has diagonals Lower, Main, Upper
n = length(Force);
x = zeros(n, 1);

% Forward elimination
for i = 2:n
    m = Lower(i) / Main(i-1);
    Main(i) = Main(i) - m * Upper(i-1);
    Force(i) = Force(i) - m * Force(i-1);
end

% Back substitution
x(n) = Force(n) / Main(n);
for i = n-1:-1:1
    x(i) = (Force(i) - Upper(i) * x(i+1)) / Main(i);
end
end