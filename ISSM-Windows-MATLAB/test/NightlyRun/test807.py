#Test Name: SquareShelfLevelsetMeltingSSA2d
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

import numpy as np

md = triangle(model(), '../Exp/Square.exp', 50000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

x = md.mesh.x
xmin = min(x)
xmax = max(x)
Lx = (xmax - xmin)
alpha = 2. / 3.
md.mask.ice_levelset = np.float_((x - alpha * Lx) > 0) - np.float_((x - alpha * Lx) < 0)

#Do not kill ice bergs as all is floating
md.levelset.kill_icebergs = 0

md.timestepping.time_step = 10
md.timestepping.final_time = 30

#Transient
md.transient.isstressbalance = True
md.transient.ismasstransport = True
md.transient.issmb = True
md.transient.isthermal = False
md.transient.isgroundingline = False
md.transient.ismovingfront = True

md.calving.calvingrate = np.zeros((md.mesh.numberofvertices))
md.frontalforcings.meltingrate = 10000 * np.ones((md.mesh.numberofvertices))
md.levelset.spclevelset = np.nan * np.ones((md.mesh.numberofvertices))

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3']
field_tolerances = [1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11,
                    2e-11, 2e-11, 2e-11, 1e-11, 1e-11, 1e-11, 5e-11,
                    2e-11, 2e-11, 2e-11, 1e-11, 1e-11, 1e-11, 5e-11]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].MaskIceLevelset,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].MaskIceLevelset,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].MaskIceLevelset]
