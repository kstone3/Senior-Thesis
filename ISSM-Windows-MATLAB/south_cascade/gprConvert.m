% Load the exported travel time data
fileID = fopen('C:\Users\bookn\Downloads\Senior-Thesis\ISSM-Windows-MATLAB\south_cascade\Data\LINE07.1.mat', 'r');
data = textscan(fileID, '%f', 'Delimiter', '\n');  % Adjust format as needed
fclose(fileID);
 % Replace with your actual filename

% Assuming 'travel_time_ns' is the variable holding travel time data
travel_time_s = travel_time_ns * 1e-9;  % Convert nanoseconds to seconds

% Define radar velocity (use your calibrated value)
velocity_m_per_ns = 0.168;

% Calculate depth
depth_m = (velocity_m_per_ns * travel_time_s) / 2;

% Display and plot
disp('Calculated Depths (in meters):');
disp(depth_m);
plot(depth_m, 'o-');
xlabel('GPR Trace Index');
ylabel('Depth (meters)');
title('Depth Profile from GPR Data');
grid on;
