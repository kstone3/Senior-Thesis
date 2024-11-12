steps=[1:4];
addpath('../bin/')
addpath('../execution/')
if any(steps==1)
    md=triangle(model, 'Data/south_cascade_glacier.exp',100);
    %plotmodel(md, 'data', 'mesh');
    save southCascade md
end
if any(steps==2)
    md=loadmodel('southCascade');
    md=setmask(md, '','');
    save southCascade md
end
if any(steps==3)
    md=loadmodel('southCascade');
    md=parameterize(md, 'southCascade.par');
    md=extrude(md,3,2);
    plotmodel(md, 'data', md.geometry.thickness)
    save southCascade md
end
if any(steps==4)
    md=loadmodel('southCascade');
    
    % Control general
    md.inversion.iscontrol = 1;
    md.inversion.nsteps = 5; % Lower the number of steps
    md.inversion.maxiter_per_step = 3 * ones(md.inversion.nsteps, 1); % Fewer iterations per step
    md.inversion.step_threshold = 0.99 * ones(md.inversion.nsteps, 1);
    md.verbose = verbose('solution', true, 'control', true);
    md.settings.solver_residue_threshold = 1e-05;  % or even 1e-04 for testing

    % Cost functions - experiment with additional cost functions for variety
    md.inversion.cost_functions = [101 103];
    md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices, 2);
    md.inversion.cost_functions_coefficients(:, 2) = 1;
    md.inversion.cost_functions_coefficients(:, 1) = 20;


    % Controls - adjust friction range
    md.inversion.control_parameters = {'FrictionCoefficient'};
    md.inversion.gradient_scaling = 1 * ones(md.inversion.nsteps, 1);
    md.inversion.min_parameters = 20 * ones(md.mesh.numberofvertices, 1);
    md.inversion.max_parameters = 80 * ones(md.mesh.numberofvertices, 1);
    md.basalforcings.groundedice_melting_rate = 1e-3; % Small baseline value
    md.basalforcings.floatingice_melting_rate = 1e-3; % Small baseline value
    md.masstransport.min_thickness=10;
    % Boundary Conditions
    % md.stressbalance.spcvx(:) = 0;
    % md.stressbalance.spcvy(:) = 0;

    % Additional parameters
    md.stressbalance.restol = 0.01;
    md.stressbalance.reltol = 0.5;
    md.stressbalance.abstol = 1e-3;

    %plotmodel(md, 'data', md.geometry.base-md.mesh.z);
    md=setflowequation(md, 'FS', 'all');
    md.cluster=generic('name',oshostname);
    md=solve(md,'Stressbalance');
    save southCascade md
end 
