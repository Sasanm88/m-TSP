import julia 
j = julia.Julia()
from julia import HGSmTSP 

import numpy as np 
from timeit import default_timer as timer



n = 500
m = 50


x = np.random.rand(n) * 1000
y = np.random.rand(n) * 1000
dist_mtx = np.sqrt((x[:,None] - x[None,:])**2 + (y[:,None] - y[None,:])**2)
dist_mtx = np.round(dist_mtx, 2)



# dist_mtx = [[   0.  ,  258.72,  823.04,  417.29,  664.02,  526.03,  644.06,
#         574.37,  616.02,  192.44],
#     [ 258.72,    0.  , 1009.86,  419.67,  758.83,  624.72,  753.46,
#         812.58,  773.15,  450.5 ],
#     [ 823.04, 1009.86,    0.  ,  696.85,  401.73,  455.61,  366.33,
#         893.78,  254.7 ,  732.85],
#     [ 417.29,  419.67,  696.85,    0.  ,  360.49,  246.91,  367.91,
#         918.38,  442.16,  525.24],
#     [ 664.02,  758.83,  401.73,  360.49,    0.  ,  138.9 ,   45.03,
#         991.52,  194.38,  674.54],
#     [ 526.03,  624.72,  455.61,  246.91,  138.9 ,    0.  ,  128.82,
#         889.03,  203.76,  548.64],
#     [ 644.06,  753.46,  366.33,  367.91,   45.03,  128.82,    0.  ,
#         952.05,  149.61,  643.74],
#     [ 574.37,  812.58,  893.78,  918.38,  991.52,  889.03,  952.05,
#         0.  ,  836.54,  401.23],
#     [ 616.02,  773.15,  254.7 ,  442.16,  194.38,  203.76,  149.61,
#         836.54,    0.  ,  572.16],
#     [ 192.44,  450.5 ,  732.85,  525.24,  674.54,  548.64,  643.74,
#         401.23,  572.16,    0.  ]]

# x = [382.56689776, 124.5875941 , 866.14750382, 227.41418686,
#     464.48309819, 433.36472507, 501.30506347, 910.30607067,
#     633.34193431, 574.89545149]
# y = [244.08118451, 224.57917736, 910.07201838, 631.45767299,
#     903.02749121, 767.65655719, 877.10289439,  17.38521931,
#     806.75055032, 237.58612663]

# dist_mtx = np.array(dist_mtx)
# x = np.array(x)
# y = np.array(y)


print(np.array2string(dist_mtx, separator=', '))
print(np.array2string(x, separator=', '))
print(np.array2string(y, separator=', '))


coordinates = np.vstack((x, y)).transpose()
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
    
