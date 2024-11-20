%estimated_base=load('Data/estimated_base.mat').estimated_base;
surface=load('Data/surface2.mat').elevationData;
x=load('Data/x2.mat').x;
y=load('Data/y2.mat').y;

velx_coords=ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','x');
vely_coords=flip(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','y'));
velx=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vx')');
vely=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vy')');
vel=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','v')');

pos12=find(isnan(velx));
velx(pos12)=0;
pos13=find(isnan(vely));
vely(pos13)=0;
pos14=find(isnan(vel));
vel(pos14)=0;

%pos18=find(abs(velx)<0.1);
%velx(pos18)=0;
%pos19=find(abs(vely)<0.1);
%vely(pos19)=0;
%pos20=find(abs(vel)<0.1);
%vel(pos20)=0;

velx = smoothdata(velx, 'movmean', 3);
vely = smoothdata(vely, 'movmean', 3);
vel = smoothdata(vel, 'movmean', 3);


md.initialization.vx=InterpFromGridToMesh(velx_coords,vely_coords,velx,md.mesh.x,md.mesh.y,0);
md.initialization.vy=InterpFromGridToMesh(velx_coords,vely_coords,vely,md.mesh.x,md.mesh.y,0);
md.initialization.vz=zeros(md.mesh.numberofvertices,1);
md.initialization.vel=InterpFromGridToMesh(velx_coords,vely_coords,vel,md.mesh.x,md.mesh.y,0);
pos20=find(md.initialization.vx<0.1 &md.initialization.vx>-0.1);
md.initialization.vx(pos20);
pos21=find(md.initialization.vy<0.1 &md.initialization.vy>-0.1);
md.initialization.vy(pos21);
pos22=find(md.initialization.vel<0.1 &md.initialization.vel>-0.1);
md.initialization.vel(pos22);

md.geometry.surface=InterpFromGridToMesh(x,y,surface,md.mesh.x,md.mesh.y,2000);

%grad_x = gradient(md.geometry.surface, 1);  % Approximate gradient in x-direction
%grad_y = gradient(md.geometry.surface')';
%slope_magnitude = sqrt(grad_x.^2 + grad_y.^2);
%thickness = (md.initialization.vel ./ (1e-24 * (917 * 9.81 * slope_magnitude).^3)).^(1/3);

[bed_elev, R] = geotiffread('Data/kriging_bed_topo_large.tif');
bed_elev=double(bed_elev);
x = R.XWorldLimits(1) + (0:size(bed_elev, 2) - 1) * R.CellExtentInWorldX;
y = flip(R.YWorldLimits(2) - (0:size(bed_elev, 1) - 1) * R.CellExtentInWorldY);
bed_elev=flipud(bed_elev);

md.geometry.base=InterpFromGridToMesh(x,y,bed_elev,md.mesh.x,md.mesh.y,2000);

pos9=find(isnan(md.geometry.base));
md.geometry.base(pos9)=0;
pos10=find(isinf(md.geometry.base));
md.geometry.base(pos10)=3000;
pos14=find(md.geometry.base>md.geometry.surface);
md.geometry.base(pos14)=md.geometry.surface(pos14);

md.geometry.thickness=md.geometry.surface-md.geometry.base;
%pos15=find(md.geometry.thickness>300);
%md.geometry.thickness(pos15)=300;

% Set minimum thickness to avoid instability
%pos_invalid = find(md.geometry.thickness <= 0.4);
%md.geometry.thickness(pos_invalid) = 0; % Set a small positive value
%md.geometry.base(pos_invalid) = md.geometry.surface(pos_invalid) - md.geometry.thickness(pos_invalid);

pos_invalid = find(md.geometry.thickness <= 1);
md.mask.ice_levelset = ones(md.mesh.numberofvertices, 1); % Default ice-covered
% Replace zero thickness with a small positive value
md.geometry.thickness(pos_invalid) = 0; % Small positive value
md.geometry.base= md.geometry.surface- md.geometry.thickness;
md.geometry.bed=md.geometry.base;

% Update ice mask
md.mask.ice_levelset(pos_invalid) = -1; % Mark as ice-free

%md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness;
surface_temperature = 263.15; % Surface temperature in Kelvin (~-10°C)
geothermal_gradient = 0.025;  % Geothermal gradient in °C/m (25 °C/km)

% Calculate depth and initialize temperature
%depth = md.geometry.thickness; % Depth from surface
%md.initialization.temperature = surface_temperature + geothermal_gradient * depth;

% Cap temperatures for physical realism
%md.initialization.temperature = min(md.initialization.temperature, 273.15); % Cap at freezing point
%md.initialization.temperature = max(md.initialization.temperature, 223.15); % Min at ~-50°C

md.initialization.temperature=269.65*ones(md.mesh.numberofvertices,1);

md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
%md.materials.rheology_B=paterson(md.initialization.temperature);

%TEMPORARY MASS BALANCE
md.smb.mass_balance = ones(md.mesh.numberofvertices, 1);
%md.smb.mass_balance(pos_invalid)=0;

md.friction.coefficient = ones(md.mesh.numberofvertices, 1) * 2;
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);
md.friction.coefficient(pos_invalid) = 1;

md.inversion=m1qn3inversion();
md.inversion.vx_obs=md.initialization.vx;
md.inversion.vy_obs=md.initialization.vy;
md.inversion.vel_obs=md.initialization.vel;
md.inversion.vz_obs=md.initialization.vz;

md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
md.masstransport.spcthickness = NaN*ones(md.mesh.numberofvertices,1);

%NEED TO SET ICE BOUNDARY TO 0, NEED TO CALCULATE ICE BOUNDARY
md.mask.ocean_levelset = ones(md.mesh.numberofvertices, 1);
%md.mask.ice_levelset = ones(md.mesh.numberofvertices, 1); % Default ice-covered
%no_ice_vertices = find(md.geometry.thickness <= 0); % Points with zero thickness
%md.mask.ice_levelset(no_ice_vertices) = -1; % Set these points as ice-free

md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);

md.stressbalance.spcvx = nan(md.mesh.numberofvertices,1); % nan = no constraint
md.stressbalance.spcvy = nan(md.mesh.numberofvertices,1);
md.stressbalance.spcvz = nan(md.mesh.numberofvertices,1);
md.stressbalance.referential = nan(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices, 3);
md=SetIceSheetBC(md);