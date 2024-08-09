
clc
clear
close
dt = 5;             % Time step [s]
t_end = 240;       % End time [s]
t = 0:dt:t_end;     % Discrete times
TS = length(t);     % Number of time steps

% Define the set of elements
H_range = [7957, 15914, 23871, 31828, 39785];

% Generate random indices
indices = randi(length(H_range), 1, TS);

% Create the random array
random_array = H_range(indices);
disp(random_array);
for i = 1: length(random_array)
    % Calculate the starting index for the current iteration
    startIndex = (i - 1) * dt + 1
    
    % Calculate the ending index for the current iteration
    endIndex = min(startIndex + dt - 1, t_end)
    HVector(startIndex:endIndex) = random_array(i);
    
end

% Display the resulting vector
plot(HVector)
xlim([0,t_end])
xticks([0:t_end/4:t_end]);
xlabel('Time [s]')
ylabel('H Amplitude, [A/m]');
title('Magnetic Field');