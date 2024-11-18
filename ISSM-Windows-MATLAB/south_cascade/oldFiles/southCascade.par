%Code to load in and calculate elevationData
%[elevationData, R] = geotiffread('Data/elevations1970_crs3413.tif');
%elevationData=double(elevationData);
%x = R.XWorldLimits(1) + (0:size(elevationData, 2) - 1) * R.CellExtentInWorldX;
%y = R.YWorldLimits(2) - (0:size(elevationData, 1) - 1) * R.CellExtentInWorldY;
%y = flip(y);               % Flip y-coordinates to be in increasing order
%writematrix(y,'Data/y.csv');
%writematrix(x,'Data/x.csv');
%save('Data/x.mat','x');
%save('Data/y.mat','y');
%elevationData = flipud(elevationData); %also need to flip elevation data 
%elevationData(elevationData == 0) = NaN;
%save('Data/surface.mat', 'elevationData');

%Code to calculate estimated ice thickness based on slope from elevationData
%dx = R.CellExtentInWorldX; % x spacing from georeferencing info
%dy = R.CellExtentInWorldY; % y spacing from georeferencing info
%[dz_dx, dz_dy] = gradient(elevationData, dx, dy);
%slope = atan(sqrt(dz_dx.^2 + dz_dy.^2));
%estimated_thickness=100000 ./ (917*9.81*sin(slope));
%pos=find(estimated_thickness<0);
%estimated_thickness(pos)=0;
%pos1=find(estimated_thickness>200);
%estimated_thickness(pos1)=200;
%disp(pos);
%disp(['Minimum thickness: ', num2str(min(estimated_thickness))]);
%disp( all(isnan(estimated_thickness(:))));
%if ~(all(estimated_thickness(:) == 0) || all(isnan(estimated_thickness(:))))
%    disp('Saving');
%    estimated_base=elevationData-estimated_thickness;
%    writematrix(estimated_base, 'Data/estimated_base.csv');
%    save('Data/estimated_base.mat','estimated_base');
%end

%Loads in estimated_base and x and y coords from code above
estimated_base=load('Data/estimated_base.mat').estimated_base;
surface=load('Data/surface.mat').elevationData;
x=load('Data/x.mat').x;
y=load('Data/y.mat').y;

%Loads in velocity
velx_coords=ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','x');
vely_coords=flip(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','y'));
velx=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vx')');
vely=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','vy')');
vel=flipud(ncread('../../Data/SouthCascadeData/ALA_G0120_1985.nc','v')');

%pos2=find(surface<estimated_base);
%surface(pos2)=estimated_base(pos2)+1;

pos3=find(isnan(surface));
surface(pos3)=0;
pos4=find(isnan(estimated_base));
estimated_base(pos4)=0;

pos5=find(isnan(velx));
velx(pos5)=0;
pos6=find(isnan(vely));
vely(pos6)=0;
pos7=find(isnan(vel));
vel(pos7)=0;

%writematrix(velx_coords, 'Data/velx.csv');
%writematrix(vely_coords, 'Data/vely.csv');

%disp(size(velx_coords));
%disp(size(vely_coords));
%disp(size(velx));
%disp(size(vely));
%disp(size(vel));

%Set model parameters
%md.miscellaneous.name='South Cascade';
%md.geometry.base=InterpFromGridToMesh(x,y,estimated_base,md.mesh.x,md.mesh.y,0);
md.geometry.surface=InterpFromGridToMesh(x,y,surface,md.mesh.x,md.mesh.y,0);
%md.geometry.thickness=md.geometry.surface-md.geometry.base;
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

%Set velocity parameters
md.inversion.vx_obs  = InterpFromGridToMesh(velx_coords,vely_coords,velx,md.mesh.x,md.mesh.y,0);
md.inversion.vy_obs  = InterpFromGridToMesh(velx_coords,vely_coords,vely,md.mesh.x,md.mesh.y,0);
md.inversion.vel_obs = InterpFromGridToMesh(velx_coords,vely_coords,vel,md.mesh.x,md.mesh.y,0);
md.initialization.vx = md.inversion.vx_obs;
md.initialization.vy = md.inversion.vy_obs;
md.initialization.vel= md.inversion.vel_obs;

grad_x = gradient(md.geometry.surface, 1);  % Approximate gradient in x-direction
grad_y = gradient(md.geometry.surface')';
slope_magnitude = sqrt(grad_x.^2 + grad_y.^2);
thickness = (md.inversion.vel_obs ./ (1e-24 * (917 * 9.81 * slope_magnitude).^3)).^(1/3);
pos8=find(thickness<0);
thickness(pos8)=0;
pos9=find(isnan(thickness));
thickness(pos9)=0;
pos10=find(isinf(thickness));
thickness(pos10)=200;
pos11=find(thickness>200);
thickness(pos11)=200;
md.geometry.base=md.geometry.surface-thickness;
md.geometry.thickness=thickness;

%Friction parameters
md.friction.coefficient=sqrt(1e+10/(917*9.81))*ones(md.mesh.numberofvertices,1).*sqrt(1/421.75);%sqrt(1./md.geometry.thickness);
md.friction.p=ones(md.mesh.numberofelements,1); % exp of tau_b
md.friction.q=ones(md.mesh.numberofelements,1); %(neg) exp on N

%Boundary Conditions
md.stressbalance.spcvx = nan(md.mesh.numberofvertices,1); % nan = no constraint
md.stressbalance.spcvy = nan(md.mesh.numberofvertices,1);
md.stressbalance.spcvz = nan(md.mesh.numberofvertices,1);
md.stressbalance.referential = nan(md.mesh.numberofvertices,6);