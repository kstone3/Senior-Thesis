%Test Name: ThermalGeothermalFlux
% This file can be run to check that the geothermal flux in simple conduction is correctly modeled.
% There is no velocity (no advection) the only thermal boundary conditions are an imposed temperature
% at upper surface and an impose flux at its base. The result must be a linear temperature from the upper to the lower
% surface with an imposed slope (Geothermal flux). if it is not the case, something is thermal modeling has been changed...
printingflag=false;

md=model();
md=triangle(md,'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareThermal.par');
md=extrude(md,11,1.);
md=setflowequation(md,'HO','all');

pos2=find(md.mesh.elementonsurface); md.thermal.spctemperature(md.mesh.elements(pos2,4:6))=0.;
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux(:)=0.1; %100mW/m^2

%analytical results
%the result is linear with depth and is equal to 0 on the upper surface (See BC)
%d2T/dz2=0  -k*dT/dz(bed)=G  T(surface)=0  => T=-G/k*(z-surface)
md.initialization.temperature=-0.1/md.materials.thermalconductivity*(md.mesh.z-md.geometry.surface); %G=0.1 W/m2

%modeled results
md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Thermal');

%plot results
comp_temp=md.results.ThermalSolution.Temperature;
relative=abs((comp_temp-md.initialization.temperature)./md.initialization.temperature)*100.;
relative(find(comp_temp==md.initialization.temperature))=0.;
plotmodel(md,'data',comp_temp,'title','Modeled temperature [K]','data',md.initialization.temperature,'view',3,...
	'title','Analytical temperature','view',3,'data',comp_temp-md.initialization.temperature,...
	'title','Absolute error [K]','view',3,'data',relative,'title','Relative error [%]','view',3,...
	'figposition','mathieu','FontSize#all',20)
if printingflag,
	set(gcf,'Color','w')
	printmodel('thermalgeothermalflux','png','margin','on','marginsize',25,'frame','off','resolution',0.7,'hardcopy','off');
	system(['mv thermalgeothermalflux.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Thermal ']);
end

%Fields and tolerances to track changes
field_names     ={'GeothermalFluxTemperature'};
field_tolerances={1e-13};
field_values    ={comp_temp};
