tic
clear all
clc
close all
%% Geometry Creation
vol_t = 126; % [mm3] Tumor Volume
% Conversion Factor for mm3 to m3 is 1e-9
vol_t_m3 = vol_t*1e-9; % [m3] Tumor Volume
RT = ((3*vol_t_m3)/(4*pi))^(1/3); % [m] Radius  of the Tumor
RH = 3*RT; % [m] Radius of the Surrounding Healthy Tissue
ST = 0.0005; % [m] Skin Thickness
RS = RH+ST; % [m] Radius of Skin Domain
R1 = [3;4;0;1;1;0;-1;-1;1;1]; % A big rectangle used to cut the geometry
C1 = [1;0;RH-RT;RT]; % Tumor Domain just below the Skin Layer
C2 = [1;0;0;RH]; % Muscle Domain
C3 = [1;0;0;RS]; % Skin Domain
C1 = [C1;zeros(length(R1) - length(C1),1)]; % Adding additional zeros for matrix concatenation
C2 = [C2;zeros(length(R1) - length(C2),1)]; % Adding additional zeros for matrix concatenation
C3 = [C3;zeros(length(R1) - length(C3),1)]; % Adding additional zeros for matrix concatenation
gd = [R1, C1, C2, C3]; % Shape combining
ns = char('R1','C1', 'C2', 'C3')'; % Name creationg for each shape
sf = '(C1+C2+C3)*R1'; % Set formula for final geometric shape

Domain_edges = decsg(gd,sf,ns);
Domain = geometryFromEdges(Domain_edges);
CPT_DM = createpde; % A dummy pde to store geometry and mesh
CPT_DM.Geometry = Domain;
generateMesh(CPT_DM,"Hmax",3E-04,"GeometricOrder","linear");
figure(1)
subplot(2,3,1)
pdegplot(CPT_DM,EdgeLabels="on",FaceLabels="on")
axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
title("Geometric Domain-Pre")
subplot(2,3,4)
pdemesh(CPT_DM); 
axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
title("RE Mesh-Pre")

%% Model Parameters

%%% Ambient conditions
htc = 10;  % [W/m2/K] Heat Transfer Coefficient
T_amb = 24; % [Degree C] Ambient Temperature

%%% Blood Properties
rho_b = 1050; % [kg/m3] Blood Density
cp_b = 3617; % [J/kg/K] Blood Specific Heat
k_b = 0.52; % [W/m/K] Blood Thermal Conductivity
T_b = 37; % [Degree C] Blood Arterial Temperature
T_b_K = T_b + 273.15; % [K] Blood Arterial Temperature

%%% Skin Properties
rho_s = 1109; % [kg/m3] Skin Density
cp_s = 3391; % [J/kg/K] Skin Specific Heat
k_s = 0.37; % [W/m/K] Skin Thermal Conducitivity
w_s = 1.96e-3; % [1/s] Skin Blood Perfusion Rate
Q_s = 1830; % [W/m3] Skin Metabolic Heat Generation Rate

%%% Muscle Properties
rho_m = 1090; % [kg/m3] Muscle Density
cp_m = 3421; % [J/kg/K] Muscle Specific Heat
k_m = 0.49; % [W/m/K] Muscle Thermal Conducitivity
w_m = 6.72e-4; % [1/s] Muscle Blood Perfusion Rate
Q_m = 992; % [W/m3] Muscle Metabolic Heat Generation Rate

%%% Tumor Properties
rho_t = 1050; % [kg/m3] Tumor Density
cp_t = 3500; % [J/kg/K] Tumor Specific Heat
k_t = 0.50; % [W/m/K] Tumor Thermal Conductivity
w_t = 8.33e-4; % [1/s] Tumor Blood Perfusion Rate
Q_t = 5790; % [W/m3] Tumor Metabolic Heat


%%% Magnetic Field Parameters
H = 10000; % [A/m] Magnetic Field Amplitude H = 10 kA/m
f = 1e5; % [Hz] Magnetic Field Frequency f = 100 kHz

%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 1.6e-8; % [m] MNP Diameter d = 16 nm
delta_l = 1e-9; % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l; % [m] MNP Hydrodynamic Diameter
rho_mnp = 5180; % [kg/m3] MNP Density
cp_mnp = 670; % [J/kg/K] MNP Specific Heat
Md = 4.46e5; % [A/m] Domain Magnetization Value Md = 446 kA/m
Keff = 20000; % [J/m3] Effective Anisotropy Constant Keff = 20 kJ/m3
MW_mnp = 0.231533; % [kg/mol] Molecular Weight of MNP MW = 231.533 g/mol

%%% MNP Dose and Magnetic Fluid Parameters
MNP_dose = 4; % [kg/m3] MNP Dose MNP mass used per unit vol of tumor = 4 mg/cm3
% Conversion 1 mg/cm3 = 1 kg/m3
MNP_mass = MNP_dose*vol_t_m3; % [kg] MNP mass injected into the tumor; Varies with MNP Dose
MNP_vol = MNP_mass/rho_mnp; % [m3] volume of MNP in the Magnetic Fluid
MF_inj_vol = 1; % [ml, milliliters] Magnetic Fluid Volume injected into the Tumor
% Converion Factor for ml to m3 is 1e-6 
MF_vol = MF_inj_vol*1e-6; % [m3] Magnetic Fluid Volume Injected in to the Tumor
rho_cf = 1000; % [kg/m3] Density of the Carrier Fluid
mu_cf = 1e-3; % [Pa.s] Viscocity of the Carrier Fluid
MNP_vol_frac = MNP_vol/MF_vol; %[-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7; % [H/m] Permeability of Free Space
kB = 1.38e-23; % [J/K] Blotzmann Constant
MNP_svol = pi*(d_mnp^3)/6; % [m3] Solid Volume of Individual Mangetic Nanoparticle
MNP_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9; % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f; % [rad/s] Angular frequency of applied magnetic field
alpha = 0.55; % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Keff*MNP_svol/(kB*T_b_K); % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*MNP_hvol)/(kB*T_b_K); % [s] Brownian Relaxation Time of MNP
tauN = sqrt(pi)*tau0*exp(gamma)/(2*sqrt(gamma)); % [s] Neel Relaxation Time of MNP
tauE = tauB*tauN/(tauB+tauN); % [s] Effective Relaxation Time of MNP
zeta = mu0*Md*H*MNP_svol/(kB*T_b_K); % [-] Langevin parameter
X_i = mu0*(Md^2)*MNP_vol_frac*MNP_svol/(3*kB*T_b_K); % [-] Initial Magnetic Susceptibility
X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta; % [-] Equilibrium Magnetic Susceptibility
X_L = X_0*omega*tauE/(1+(omega*tauE)^2); % [-] Loss Component of Magnetic Susceptibility
P = pi*mu0*f*X_L*H^2; % [W/m3] Heat Generation rate of MNPs in MF
SAR = P/(rho_mnp*MNP_vol_frac); % [W/kg] MNP Specific Absorption Rate
SAR_grams = SAR/1000; % [W/g] MNP Specific Absorption Rate

%%% Heat Source by MNP in Tissue

MNP_conc = MNP_mass/vol_t_m3; % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
Q_MNP = alpha*MNP_conc*SAR; % [W/m3] Heat Generation by MNP in Tissue


%%

modelSS = createpde(1);
modelSS.Geometry = CPT_DM.Geometry;
modelSS.Mesh = CPT_DM.Mesh;
subplot(2,3,2)
pdegplot(modelSS.Geometry,EdgeLabels="on",FaceLabels="on")
axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
title("Geometric Domain-SS")
subplot(2,3,5)
pdemesh(modelSS); 
axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
title("RE Mesh-SS")

specifyCoefficients(modelSS, "m",0,"d",0,"c",k_m,"a",rho_b*cp_b*w_m,"f",rho_b*cp_b*w_m*T_b+Q_m, 'Face',1);
specifyCoefficients(modelSS, "m",0,"d",0,"c",k_s,"a",rho_b*cp_b*w_s,"f",rho_b*cp_b*w_s*T_b+Q_s, 'Face',2);
specifyCoefficients(modelSS, "m",0,"d",0,"c",k_t,"a",rho_b*cp_b*w_t,"f",rho_b*cp_b*w_t*T_b+Q_t+Q_MNP, 'Face',3);
applyBoundaryCondition(modelSS,"dirichlet","Edge",9:10,"h",1,"r",T_amb);
applyBoundaryCondition(modelSS,"neumann","Edge",1:4,"q",0,"g", 0);

% thermalBC(modelSS,"Edge",[10,11],"Temperature",T_b);
% thermalBC(modelSS,"Edge",[9,10],"ConvectionCoefficient",htc,"AmbientTemperature",T_amb);
SST = solvepde(modelSS);

figure 
pdeplot(modelSS,"XYData",SST.NodalSolution,"Contour","on")
axis equal
title('Steady-State Temperature')
% %% Transient Simulation
% 
% modelTS = createpde("thermal","transient-axisymmetric");
% modelTS.Geometry = CPT_DM.Geometry;
% modelTS.Mesh = CPT_DM.Mesh;
% subplot(2,3,3)
% pdegplot(modelTS.Geometry,EdgeLabels="on",FaceLabels="on")
% axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
% title("Geometric Domain-TS")
% subplot(2,3,6)
% pdemesh(modelTS); 
% axis([-.25*RH, 1.25*RH, -1.25*RH, 1.25*RH]);
% title("RE Mesh-TS")
% thermalProperties(modelTS,'Face',1, "ThermalConductivity",k_m, ...
%                                     "MassDensity",rho_m, ...
%                                     "SpecificHeat",cp_m);
% thermalProperties(modelTS,'Face',2, "ThermalConductivity",k_t, ...
%                                     "MassDensity",rho_t, ...
%                                     "SpecificHeat",cp_t);
% internalHeatSource(modelTS, @(location,state) rho_b*cp_b*w_m*(T_b- state.u), 'Face',1);
% internalHeatSource(modelTS, @(location,state) rho_b*cp_b*w_t*(T_b- state.u), 'Face',2);
% 
% internalHeatSource(modelTS, Q_m, 'Face',1);
% internalHeatSource(modelTS, Q_t, 'Face',2);
% internalHeatSource(modelTS, Q_MNP, 'Face',2);
% 
% thermalBC(modelTS,"Edge",[6,7],"Temperature",T_b);
% % IC = @(location,state) SST.Temperature;
% thermalIC(modelTS, SST.Temperature);
% 
% tlist = 0:5:1200;
% 
% TST = solve(modelTS,tlist);
% 
% figure 
% pdeplot(modelTS,"XYData",TST.Temperature(:,end),"Contour","on")
% axis equal
% title('Post-treatment Temperature')
% 
% toc