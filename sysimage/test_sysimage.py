
from julia.api import LibJulia 
api = LibJulia.load()
api.sysimage = "sysimage/sys_hgsmtsp.so"  # path to sys_hgsmtsp.so
api.init_julia()

from julia import HGSmTSP

import numpy as np 
x = np.random.rand(100) * 1000
y = np.random.rand(100) * 1000
dist_mtx = np.sqrt((x[:,None] - x[None,:])**2 + (y[:,None] - y[None,:])**2)


coordinates = np.concatenate((x, y)).reshape(-1, 2)
HGSmTSP.solve_mTSP(3, dist_mtx, coordinates)
