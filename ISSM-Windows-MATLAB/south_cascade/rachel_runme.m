%% Add ISSM commands to current path
addpath /Users/rachelmiddleton/Workspace/ISSM/ISSM-macOS-Intel-MATLAB/bin /Users/rachelmiddleton/Workspace/ISSM/ISSM-macOS-Intel-MATLAB/lib /Users/rachelmiddleton/Workspace/ISSM/ISSM-macOS-Intel-MATLAB/share/proj

steps=[];
if any(steps==0) 
    %% Negribreen Geometry - create an .exp file from geometry
    
    % the contour file / outline of glacier
    contour_file = 'Geometry/contour_negri_glims_less_UTM_20170707_v2.dat';
    
    % Create the outline file for meshing
    outline = load(contour_file);
    
    
    
    % write NegribreenOuline.exp
    filename = 'Geometry/NegribreenOutline.exp';
    fid = fopen(filename, 'wt');
    
    fprintf(fid,'## Name:NegribreenOutline\n');
    fprintf(fid,'## Icon:0\n');
    fprintf(fid,'# Points Count Value\n');
    fprintf(fid,'%d 1.000000\n',length(outline));
    fprintf(fid,'# X pos Y pos\n');
    
    for i = 1:length(outline)
        fprintf(fid, '%.6f %.6f\n', outline(i,1), outline(i,2));
    end
    
    fclose(fid);
    
    % Create the mesh
    md=model;
    md=triangle(md,'Geometry/NegribreenOutline.exp',400);
    plotmodel(md,'data','mesh');
    
    save NegribreenModel md
end


if any(steps==1) 
	disp('   Step 1: Read mesh');
	loadmodel NegribreenModel;

    

	%Get observed velocity field on mesh nodes

    % Use the velocity data located here: 
    % 'Velocity/subset_1_of_S1A_IW_GRDH_1SSH_20150802T154516_20150814T154541_007086_009A5A_DFBB_Orb_Stack_vel_corr1_hfr4_rww256_data/vector_data/Velocity.csv'
    vel_data_path = 'Velocity/subset_1_of_S1A_IW_GRDH_1SSH_20150802T154516_20150814T154541_007086_009A5A_DFBB_Orb_Stack_vel_corr1_hfr4_rww256_data/vector_data/Velocity.csv'
    vel_data = readtable(vel_data_path);
    lat_slv = vel_data.slv_lat_Double;
    lon_slv = vel_data.slv_lon_Double;
    vel = vel_data.velocity_Double; % this is m/d
    vel = vel.*365.242; % convert to m/y

    [x1,y1] = ll2utm(lat_slv,lon_slv,33); % hope everything else is in utm 33
    % put vel on a grid
    x_grid = unique(x1);
    y_grid = unique(y1);
    [Xq, Yq] = meshgrid(x_grid,y_grid);
    Vq = griddata(x1, y1, vel, Xq, Yq, 'linear');

    vel = InterpFromGridToMesh(Xq,Yq,Vq,md.mesh.x,md.mesh.y,0); % need to figure out how to get unit vectors (vx and vy)

	% x1		= ncread(ncdata,'x1');
	% y1		= ncread(ncdata,'y1');
	% velx	= ncread(ncdata,'surfvelx');
	% vely	= ncread(ncdata,'surfvely');
	% vx		= InterpFromGridToMesh(x1,y1,velx',md.mesh.x,md.mesh.y,0);
	% vy		= InterpFromGridToMesh(x1,y1,vely',md.mesh.x,md.mesh.y,0);
	% vel		= sqrt(vx.^2+vy.^2);

	%refine mesh using surface velocities as metric
	md=bamg(md,'hmin',1200,'hmax',15000,'field',vel,'err',5);
	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,33,71);
    

	save NegribreenModel md
end 

if any(steps==2) 
	disp('   Step 2: Parameterization');
	md=loadmodel('JksMesh');

	md=setmask(md,'','');
	md=parameterize(md,'Jks.par');

	save JksPar md
end 

if any(steps==3) 
	disp('   Step 3: Control method friction');
	md=loadmodel('JksPar');

	md=setflowequation(md,'SSA','all');

	%Control general
	md.inversion.iscontrol=1;
	md.inversion.nsteps=20;
	md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
	md.inversion.maxiter_per_step=5*ones(md.inversion.nsteps,1);
	md.verbose=verbose('solution',true,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=40;
	md.inversion.cost_functions_coefficients(:,2)=1;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.gradient_scaling(1:md.inversion.nsteps)=30;
	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=generic('name',oshostname,'np',4);
	md=solve(md,'Stressbalance');

	save JksControl md
end 

if any(steps==4) 
	disp('   Plotting')
	md=loadmodel('JksControl');

	plotmodel(md,'unit#all','km','axis#all','equal',...
		'data',md.inversion.vel_obs,'title','Observed velocity',...
		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
		'colorbar#1','off','colorbartitle#2','(m/yr)',...
		'caxis#1-2',[0,7000],...
		'data',md.geometry.base,'title','Base elevation',...
		'data',md.results.StressbalanceSolution.FrictionCoefficient,...
		'title','Friction Coefficient',...
		'colorbartitle#3','(m)');
end 
