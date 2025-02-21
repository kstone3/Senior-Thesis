steps=[3];
addpath('../bin/')
addpath('../execution/')

if any(steps==1)
    domain =['Data/south_cascade_glacier.exp'];
	hinit=1000;	% element size for the initial mesh
	hmax=1000;		% maximum element size of the final mesh
	hmin=100;		% minimum element size of the final mesh
	gradation=1.3;	% maximum size ratio between two neighboring elements
	err=5;			% maximum error between interpolated and control field

	% Generate an initial uniform mesh (resolution = hinit m)
	md=bamg(model,'domain',domain,'hmax',hinit);

    velx_coords=ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','x');
    vely_coords=flip(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','y'));
    vel=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','v')');
	vel_obs=InterpFromGridToMesh(velx_coords,vely_coords,vel,md.mesh.x,md.mesh.y,0);
	pos14=find(isnan(vel_obs));
	vel_obs(pos14)=0;
	pos20=find(abs(vel)<0.1);
	vel(pos20)=0;

    md=bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',vel_obs,'err',err);
	%md=triangle(model, 'Data/south_cascade_glacier.exp',50);
	%ploting
	% plotmodel(md,'data','mesh')
	save Models/southCascadeMesh md;
end
if any(steps==2) %Parameterization #3 
	md=loadmodel('Models/southCascadeMesh');
	md = parameterize(md,'southCascade2.par');
	disp('Checking for nan')
	disp(['Surface nan: ', num2str(any(isnan(md.geometry.surface)))])
	disp(['Base nan: ', num2str(any(isnan(md.geometry.base)))])
	disp(['Thickness<0: ', num2str(any(md.geometry.thickness < 0))]) % Negative or zero thickness
	disp(['vx nan: ', num2str(any(isnan(md.initialization.vx)))])
	disp(['vy nan: ', num2str(any(isnan(md.initialization.vy)))])
	disp(['vel nan: ', num2str(any(isnan(md.initialization.vel)))])
	disp(['friction nan: ', num2str(any(isnan(md.friction.coefficient)))])
	disp(['friction<=0: ', num2str(any(md.friction.coefficient <= 0))]) % Invalid friction coefficients
	disp(['SMB Min: ', num2str(min(md.smb.mass_balance))]);
	disp(['SMB Max: ', num2str(max(md.smb.mass_balance))]);
	disp(['Min temperature: ', num2str(min(md.initialization.temperature))]);
	disp(['Max temperature: ', num2str(max(md.initialization.temperature))]);
	min_thickness = min(md.geometry.thickness);
	disp(['Min thickness: ', num2str(min(md.geometry.thickness))]);
	disp(['Number of invalid thickness values: ', num2str(sum(md.geometry.thickness < 0))]);
	disp(['Number of base > surface: ', num2str(sum(md.geometry.base > md.geometry.surface))]);
	min_friction = min(md.friction.coefficient);
	disp(['Minimum friction: ', num2str(min_friction)]);
	%md.thermal.isthermal = 0;
	% plotmodel(md, 'data', md.geometry.thickness);
	save Models/southCascadePar md;
end
if any(steps==3)
	md=loadmodel('Models/southCascadePar');
	md=extrude(md,3,1);
	%plotmodel(md, 'data', md.geometry.thickness, 'title', 'Ice Thickness');
	%plotmodel(md, 'data', md.geometry.base, 'title', 'Bed Elevation');
	%plotmodel(md, 'data', md.friction.coefficient, 'title', 'Friction Coefficient');
	%plotmodel(md, 'data', md.initialization.vx, 'title', 'Initial Velocity X');
	%plotmodel(md, 'data', md.initialization.vy, 'title', 'Initial Velocity Y');
	%plotmodel(md, 'data', md.smb.mass_balance, 'title', 'Surface Mass Balance');
	% Extract mesh data
	elements = md.mesh.elements; % Triangular elements (connectivity)
	x = md.mesh.x; % Node x-coordinates
	y = md.mesh.y; % Node y-coordinates

	% Compute edge lengths for each element
	edge1 = sqrt((x(elements(:,2)) - x(elements(:,1))).^2 + (y(elements(:,2)) - y(elements(:,1))).^2);
	edge2 = sqrt((x(elements(:,3)) - x(elements(:,2))).^2 + (y(elements(:,3)) - y(elements(:,2))).^2);
	edge3 = sqrt((x(elements(:,1)) - x(elements(:,3))).^2 + (y(elements(:,1)) - y(elements(:,3))).^2);

	% Find the longest and shortest edge for each element
	longest_edge = max([edge1, edge2, edge3], [], 2);
	shortest_edge = min([edge1, edge2, edge3], [], 2);

	% Calculate aspect ratio for each element
	aspect_ratios = longest_edge ./ shortest_edge;

	% Display statistics
	disp(['Max aspect ratio: ', num2str(max(aspect_ratios))]);
	disp(['Average aspect ratio: ', num2str(mean(aspect_ratios))]);
	% Compute areas of triangular elements
	elements = md.mesh.elements;
	x = md.mesh.x;
	y = md.mesh.y;
	area = 0.5 * abs( ...
		x(elements(:,1)) .* (y(elements(:,2)) - y(elements(:,3))) + ...
		x(elements(:,2)) .* (y(elements(:,3)) - y(elements(:,1))) + ...
		x(elements(:,3)) .* (y(elements(:,1)) - y(elements(:,2))) );

	% Identify elements with very small or negative areas
	invalid_elements = find(area <= 0);
	disp(['Number of invalid elements: ', num2str(length(invalid_elements))]);

	% Visualize problematic elements
	if ~isempty(invalid_elements)
		plotmodel(md, 'data', invalid_elements, 'title', 'Invalid Elements');
	end

	% plotmodel(md, 'data', 'mesh');
	% % Solve

	md = setflowequation(md,'SIA','all');
	md.cluster=generic('name',oshostname,'np',2);
	md.toolkits=toolkits;
	md.verbose = verbose('solution', true);
	md.stressbalance.reltol = 1e-3;  % Increase solver tolerance
	md.stressbalance.maxiter = 500;  % Allow more iterations

	% md.verbose.solution=1;
	md.timestepping.time_step = 0.1;    % Define time step (e.g., 0.1 years)
	%md.timestepping.final_time = 10;   % Define final simulation time (e.g., 10 years)
	md = solve(md, 'transient');       % Solve the transient problem
	% md.timestepping.time_step = 0;
	% md.timestepping.final_time = 0;
	% md = solve(md, 'SteadyState');	

	save Models/southCascadeSolved md;
end
if any(steps==4)
	md=loadmodel('Models/southCascadeSolved.mat');
	fieldnames(md.results.SteadystateSolution)
	times = [md.results.SteadystateSolution.time]; % Extract time points
	n_steps = length(md.results.SteadystateSolution); % Number of time steps
	filename = 'VelocityChangeOverTime.gif'; % Name of the GIF file
	% Check one time step (e.g., first time step)
	any(md.results.SteadystateSolution(1).Vel > 0) % Should return true if velocity is non-zero
	% Check across all time steps
	for t = 1:length(md.results.SteadystateSolution)
		disp(['Time step ', num2str(t), ': Max Velocity = ', num2str(max(md.results.SteadystateSolution(t).Vel))]);
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