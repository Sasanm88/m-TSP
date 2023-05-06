import julia 
j = julia.Julia()
from julia import HGSmTSP 

import numpy as np 
from timeit import default_timer as timer



n = 10
m = 3


x = np.random.rand(n) * 1000
y = np.random.rand(n) * 1000
dist_mtx = np.sqrt((x[:,None] - x[None,:])**2 + (y[:,None] - y[None,:])**2)


coordinates = np.concatenate((x, y)).reshape(-1, 2)
t0 = timer()
hgs_routes, hgs_route_lengths = HGSmTSP.solve_mTSP(m, dist_mtx, coordinates, n_iterations=1000, time_limit=2.5)
t1 = timer()
print(f"---- HGSmTSP ------")
print(f"Time: {t1-t0}")
print(f"Obj : {max(hgs_route_lengths)}")

hgs_time = t1 - t0 
# Comparison with NCE 
try:
    import nce.solver

    t0 = timer()
    routes, route_lengths = nce.solver.solve_mTSP(m, dist_mtx, num_candidates=5, perturbation=3, time_limit=hgs_time)
    t1 = timer()
    print("---- NCE ------")
    print(f"Time: {t1-t0}")
    print(f"Obj : {max(route_lengths)}")
except:
    print("NCE is not available.")
    
