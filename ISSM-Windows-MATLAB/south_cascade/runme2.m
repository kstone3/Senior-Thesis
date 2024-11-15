steps=[2:3];
addpath('../bin/')
addpath('../execution/')
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
	md=extrude(md,3,1);
	md = setflowequation(md,'FS','all');
	plotmodel(md, 'data', md.geometry.thickness);
	save Models/southCascadePar md;
end
if any(steps==3)
	md=loadmodel('Models/southCascadePar');
	% plotmodel(md, 'data', md.geometry.base);
	md.timestepping.time_step=0.01;
	md.timestepping.final_time=10;
	md.settings.solver_residue_threshold=1e-5;
	% Solve
	md.toolkits=toolkits;
	md.cluster=generic('name',oshostname,'np',2);
	md=solve(md,'Transient');

	save Models/southCascadeSolved md;
end