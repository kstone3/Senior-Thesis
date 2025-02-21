addpath('../bin/')
addpath('../execution/')
%% South Cascade Glacier Volume Change (1984-2024) using ISSM
clear; close all; clc;



%% Step 1: Load DEMs (Surface Elevation)
% disp('Loading DEMs...');
% DEM_1986 = load('Data/elevations1986_crs3413.tif'); % Load 1984 DEM
% % DEM_2024 = load('DEM_2024.asc'); % Load 2024 DEM

% % Extract x, y, and elevation
% [x_1984, y_1984, z_1984] = dem_to_xyz(DEM_1986);
% [x_2024, y_2024, z_2024] = dem_to_xyz(DEM_2024);

%% Step 2: Create ISSM Mesh
disp('Generating Mesh...');
md = triangle(model(),'Data/south_cascade_glacier.exp',500); % 500m resolution
md = setmask(md,'all',''); % Define glacier mask
disp('Setting up time stepping...');
md.timestepping.start_time = 1984;
md.timestepping.final_time = 2024; % 40-year simulation
md.timestepping.time_step = 1; % 1-year step
md = parameterize(md,'southCascade3.par');
disp('Defining flow model...');
md = setflowequation(md, 'SSA', 'all'); % Use SSA (Shallow Shelf Approximation)


%% Step 3: Assign Ice Thickness and Bed Topography
disp('Loading Ice Thickness Data...');
% thickness_1984 = load('ice_thickness_1984.txt'); % Load ice thickness
% thickness_2024 = load('ice_thickness_2024.txt'); % Load 2024 thickness (if available)

% Compute bed topography
surface=load('Data/surface2.mat').elevationData;
x=load('Data/x2.mat').x;
y=load('Data/y2.mat').y;
md.geometry.surface=InterpFromGridToMesh(x,y,surface,md.mesh.x,md.mesh.y,2000);
[bed_elev, R] = geotiffread('Data/kriging_bed_topo_large.tif');
bed_elev=double(bed_elev);
x = R.XWorldLimits(1) + (0:size(bed_elev, 2) - 1) * R.CellExtentInWorldX;
y = flip(R.YWorldLimits(2) - (0:size(bed_elev, 1) - 1) * R.CellExtentInWorldY);
bed_elev=flipud(bed_elev);

md.geometry.base=InterpFromGridToMesh(x,y,bed_elev,md.mesh.x,md.mesh.y,2000);
md.geometry.base(find(md.geometry.base>md.geometry.surface))=md.geometry.surface(find(md.geometry.base>md.geometry.surface));
md.geometry.thickness = md.geometry.surface - md.geometry.base;
% Ensure thickness is positive
% md.geometry.thickness(md.geometry.thickness <= 0) = 1; % Minimum thickness of 1m to prevent errors


%% Step 4: Apply Surface Mass Balance (SMB)
% disp('Applying SMB...');
% md.smb.mass_balance = load('smb_1984_2024.txt'); % Load SMB data (e.g., in m/yr)

%% Step 5: Set Up Time Stepping for Transient Simulation
disp('Setting up time stepping...');
md.timestepping.start_time = 1984;
md.timestepping.final_time = 2024; % 40-year simulation
md.timestepping.time_step = 1; % 1-year step

%% Step 6: Set Up Physics (Ice Flow, Mass Transport)
disp('Defining physics...');
md.transient.isstressbalance = 1;
md.transient.ismasstransport = 1;
md.transient.isthermal = 0; % No thermal evolution
md.transient.ismovingfront = 1;

disp('Checking for NaN values in model fields...');

% List of fields to check
fields_to_check = {
    md.geometry.surface, 
    md.geometry.base, 
    md.geometry.thickness, 
    md.smb.mass_balance, 
    md.materials.rheology_n, 
    md.materials.rheology_B, 
    md.friction.coefficient, 
    md.friction.p, 
    md.friction.q, 
    md.basalforcings.groundedice_melting_rate, 
    md.basalforcings.floatingice_melting_rate
};

% Iterate through fields and check for NaN values
for i = 1:length(fields_to_check)
    if any(isnan(fields_to_check{i}))
        error(['NaN values detected in model field #', num2str(i)]);
    end
end



%% Step 7: Solve the Model
disp('Running ISSM simulation...');
md = solve(md, 'Transient');

%% Step 8: Compute Volume Change
disp('Computing volume change...');
V_1984 = sum(md.results.TransientSolution(1).Thickness .* md.mesh.elements_area);
V_2024 = sum(md.results.TransientSolution(end).Thickness .* md.mesh.elements_area);
dV = V_2024 - V_1984;

disp(['Total Volume Change (1984-2024): ', num2str(dV), ' cubic meters']);

%% Step 9: Plot Results
disp('Plotting results...');
figure;
trisurf(md.mesh.elements, md.mesh.x, md.mesh.y, ...
       md.results.TransientSolution(end).Thickness - md.results.TransientSolution(1).Thickness);
title('South Cascade Glacier Volume Change (1984-2024)');
colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Thickness Change [m]');
view(3);
grid on;

disp('Simulation complete!');
