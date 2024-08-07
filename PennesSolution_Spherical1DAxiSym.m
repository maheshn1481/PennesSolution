tic
clc
clear
close all
%%% Solves the Pennes BioHeat Equation in Spherical Coordinates %%%
%%% The model is 1-Dimensional and Axisymmetric %%%
%%% The initial conditions are T(0,r) = 37 °C %%%
%%% Boundary Conditions are Symmetry at r = 0; and Fixed Temperature at r = R as T(0,R) = 37 °C %%%
%%% Uses TDMA algorithm to solve Discretized Linear Equations %%%


%%% Model Inputs
modelInput = struct;
modelInput.tumor_size = 0.005;% [m] Radius of the Tumor
modelInput.treatment_time = 1200;% End time [s]
modelInput.initial_condition = 37;% [°C] Initial Condition Temperature

%%% Magnetic Field Application status
heat_source = struct;
heat_source.type = 'Constant';% Magnetic field is on through out the treatment time
heat_source.on_time = modelInput.treatment_time;% [s] Magnetic field On time, if Magnetic field is on through out the treatment time
% heat_source.type = 'Pulsed';% Magnetic field is turned on and off
% heat_source.on_time = 10;% [s] Magnetic field On time, if Magnetic field is on through out the treatment time

%%% Blood Properties
Blood = struct;
Blood.rho = 1000;           % [kg/m3] Blood Density
Blood.cp = 3840;            % [J/kg/K] Blood Specific Heat
Blood.T = 37;               % [°C] Blood Arterial Temperature

%%% Tumor Properties
Tumor = struct;
Tumor.rho = 1045;% [kg/m3] Tumor Density
Tumor.cp = 3600;% [J/kg/K] Tumor Specific Heat
Tumor.k = 0.527;% [W/m/K] Tumor Thermal Conductivity
Tumor.w = 0.009;% [1/s] Tumor Blood Perfusion Rate
Tumor.Qm = 0;% [W/m3] Tumor Metabolic Heat

%%% Magnetic Field Parameters
Mag_Field = struct;
Mag_Field.Amplitude = 20000;% [A/m] Magnetic Field Amplitude H = 10 kA/m
Mag_Field.Frequency = 1e5;% [Hz] Magnetic Field Frequency f = 100 kHz

[Temperature,r,t,q_mnp] = pennes_temperature_solver(modelInput,heat_source,Blood,Tumor,Mag_Field);
% [Temperature,r,t,q_mnp] = pennes_temperature_solver2(modelInput,heat_source,Blood,Tumor,Mag_Field);
T_sum_time = sum(Temperature,2);      % Summation of Temperature over all Times at each r location

for i = 1:size(Temperature,2)
    GAIN(1,i) = trapz(r, 4*pi*r.^2.*Temperature(:,i));     % Integration of Temperature over the domain
end


Total_GAIN = trapz(r, 4*pi*r.^2.*T_sum_time)     % Integration of Temperature over the domain


for i = 1:size(Temperature,2)
    volume(1,i) = trapz(r, 4*pi*r.^2.*1);     % Integration of radial position over the domain (Volume of sphere)
end
toc
%%
%%% Plot the results
close all
figure('Position', [100, 100, 800, 400])
sgtitle([heat_source.type, ' Heat Source' ])
subplot(1,2,1)
plot_radius = [0,modelInput.tumor_size,5*modelInput.tumor_size];
legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
plot_rad_idx = (plot_radius/(r(2)-r(1)))+1;
plot(t,Temperature(plot_rad_idx,:), 'LineWidth',2)
xlabel('Time, t [s]');
ylabel('Temperature, T [°C]');
xlim([0,t(end)])
xticks(0:200:t(end));
legend(legend_String, Location='south')
title('Temperature Elevations');
grid on;
subplot(1,2,2)
plot_time = [0,60,120,t(end)];
legend_String = string(plot_time)+[" s"," s"," s"," s"];
plot_ind = (plot_time/(t(2)-t(1)))+1;
plot(r,Temperature(:,plot_ind),'LineWidth',2)
hold on
xlabel('Radial distance, r [m]');
xlim([0,5*modelInput.tumor_size])
xticks(0:5*modelInput.tumor_size/5:5*modelInput.tumor_size);
ylabel('Temperature, T [°C]');
legend(legend_String)
title('Spatial Temperature profile');
