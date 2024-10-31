steps=[3:4]
addpath('../bin/')
if any(steps==1)
    md=triangle(model, 'south_cascade_glacier.exp',0.1);
    %plotmodel(md, 'data', 'mesh');
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
    %disp(md.mesh.x)
    save southCascade md
end
if any(steps==4)
    md=loadmodel('southCascade');
    %disp(md.geometry.base)
    plotmodel(md, 'data', md.geometry.base);
    save southCascade md
end 
