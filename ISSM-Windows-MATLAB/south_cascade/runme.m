steps=[3]
addpath('../bin/')
if any(steps==1)
    md=triangle(model, 'south_cascade_glacier.exp',100);
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
    %disp(md.geometry.base)
    %plotmodel(md, 'data', md.geometry.base);
    %disp(md.mesh.x)
    save southCascade md
end
if any(steps==4)
    md=loadmodel('southCascade');
    %md.mesh.epsg = 32610; % Replace 32610 with your desired EPSG code
    % md = reproject(md, md.mesh.epsg);
    disp(md.mesh)
    %disp(md.geometry.base)
    %plotmodel(md, 'data', md.geometry.base);
    save southCascade md
end 
