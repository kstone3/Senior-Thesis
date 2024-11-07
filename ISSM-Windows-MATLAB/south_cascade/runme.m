steps=[4]
addpath('../bin/')
if any(steps==1)
    md=triangle(model, 'Data/south_cascade_glacier.exp',100);
    plotmodel(md, 'data', 'mesh');
    save southCascade md
end
if any(steps==2)
    md=loadmodel('southCascade');
    md=setmask(md, '','');
    %plotmodel(md, 'data', md.mask.ice_levelset)
    save southCascade md
end
if any(steps==3)
    md=loadmodel('southCascade');
    md=parameterize(md, 'southCascade.par')
    md=extrude(md,5,2)
    %disp(md.geometry.base)
    %plotmodel(md, 'data', md.inversion.vel_obs);
    %disp(md.mesh.x)
    save southCascade md
end
if any(steps==4)
    md=loadmodel('southCascade');
    md.inversion.nsteps=20;
    %plotmodel(md, 'data', md.inversion.vel_obs);
    md=setflowequation(md, 'FS', 'all')
    md=solve(md,'Masstransport');
    save southCascade md
end 
