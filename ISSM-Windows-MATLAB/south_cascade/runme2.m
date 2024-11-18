steps=[2:4];
addpath('../bin/')
addpath('../execution/')

% [elevationData, R] = geotiffread('Data/elevations1986_crs3413.tif');
% elevationData=double(elevationData);
% x = R.XWorldLimits(1) + (0:size(elevationData, 2) - 1) * R.CellExtentInWorldX;
% y = R.YWorldLimits(2) - (0:size(elevationData, 1) - 1) * R.CellExtentInWorldY;
% y = flip(y);               % Flip y-coordinates to be in increasing order
% % writematrix(y,'Data/y.csv');
% % writematrix(x,'Data/x.csv');
% save('Data/x2.mat','x');
% save('Data/y2.mat','y');
% elevationData = flipud(elevationData); %also need to flip elevation data 
% elevationData(elevationData == 0) = NaN;
% save('Data/surface2.mat', 'elevationData');

if any(steps==1)
    % domain =['Data/south_cascade_glacier.exp'];
	% hinit=10000;	% element size for the initial mesh
	% hmax=40000;		% maximum element size of the final mesh
	% hmin=5000;		% minimum element size of the final mesh
	% gradation=1.7;	% maximum size ratio between two neighboring elements
	% err=8;			% maximum error between interpolated and control field

	% % Generate an initial uniform mesh (resolution = hinit m)
	% md=bamg(model,'domain',domain,'hmax',hinit);

    % velx_coords=ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','x');
    % vely_coords=flip(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','y'));
    % velx=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vx')');
    % vely=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vy')');
    % vel=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','v')');

    % vx_obs=InterpFromGridToMesh(velx_coords,vely_coords,velx,md.mesh.x,md.mesh.y,0);
	% vy_obs=InterpFromGridToMesh(velx_coords,vely_coords,vely,md.mesh.x,md.mesh.y,0);
	% vel_obs=InterpFromGridToMesh(velx_coords,vely_coords,vel,md.mesh.x,md.mesh.y,0);

    % md=bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',vel_obs,'err',err);
	md=triangle(model, 'Data/south_cascade_glacier.exp',100);
	%ploting
	%plotmodel(md,'data','mesh')
	save Models/southCascadeMesh md;
end
if any(steps==2) %Parameterization #3 
	md=loadmodel('Models/southCascadeMesh');
	md = parameterize(md,'southCascade2.par');
	disp('Checking for nan')
	any(isnan(md.geometry.surface))
	any(isnan(md.geometry.base))
	any(md.geometry.thickness <= 0) % Negative or zero thickness
	any(isnan(md.initialization.vx))
	any(isnan(md.initialization.vy))
	any(isnan(md.initialization.vel))
	any(isnan(md.friction.coefficient))
	any(md.friction.coefficient <= 0) % Invalid friction coefficients
	max(md.initialization.vx)         % Should show a reasonable value
	max(md.initialization.vy)
	md=extrude(md,3,1);
	md = setflowequation(md,'FS','all');
	plotmodel(md, 'data', md.initialization.vel);
	save Models/southCascadePar md;
end
if any(steps==3)
	md=loadmodel('Models/southCascadePar');
	% plotmodel(md, 'data', md.geometry.base);
	% % Solve
	md.toolkits=toolkits;
	md.timestepping.time_step = 0.1;    % Define time step (e.g., 0.1 years)
	md.timestepping.final_time = 10;   % Define final simulation time (e.g., 10 years)
	md = solve(md, 'Transient');       % Solve the transient problem
	md.cluster=generic('name',oshostname,'np',2);
	% md.timestepping.time_step = 0; % Disable time stepping
	% md.timestepping.final_time = 0; % Steady state
	%md = solve(md, 'SteadyState'); % Solve in steady state

	%save Models/southCascadeSolved md;
end
if any(steps==4)
	md=loadmodel('Models/southCascadeSolved.mat')
	fieldnames(md.results.TransientSolution)
	times = [md.results.TransientSolution.time]; % Extract time points
	n_steps = length(md.results.TransientSolution); % Number of time steps
	filename = 'VelocityChangeOverTime.gif'; % Name of the GIF file
	% Check one time step (e.g., first time step)
	any(md.results.TransientSolution(1).Vel > 0) % Should return true if velocity is non-zero

	% Check across all time steps
	for t = 1:length(md.results.TransientSolution)
		disp(['Time step ', num2str(t), ': Max Velocity = ', num2str(max(md.results.TransientSolution(t).Vel))]);
	end

	% for t = 1:n_steps
	% 	figure;
	% 	plotmodel(md, 'data', md.results.TransientSolution(t).Vel, 'title', ...
	% 		['Velocity Magnitude at Time Step ', num2str(t), ' (Time = ', num2str(times(t)), ' years)']);
	% 	colorbar;

	% 	% Capture the frame and convert it to an image
	% 	frame = getframe(gcf);
	% 	img = frame2im(frame);
	% 	[imind, cm] = rgb2ind(img, 256);
		
	% 	% Write to GIF
	% 	if t == 1
	% 		imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
	% 	else
	% 		imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
	% 	end
	% 	close;
	% end
end