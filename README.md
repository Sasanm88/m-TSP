This page contains the HGA code proposed in the following paper: 

Mahmoudinazlou, S., & Kwon, C. (2023). A Hybrid Genetic Algorithm for the min-max Multiple Traveling Salesman Problem. arXiv preprint arXiv:2307.07120.

## How to Install as a Julia package?

```
cd m-TSP
julia
] add .
```

Test:
```julia
using HGSmTSP

n_vehicles = 3
n_nodes = 10
coordinates = rand(n_nodes, 2) .* 1000

dist_mtx = zeros(n_nodes, n_nodes)
for i=1:n_nodes
    for j=1:n_nodes
        dist_mtx[i, j] = sum(abs.(coordinates[i, :] - coordinates[j, :]))
    end
end

@time solve_mTSP(n_vehicles, dist_mtx, coordinates)
```

In the above `dist_mtx` is a matrix of size `n+1` by `n+1` where `n` is the number of customers. The first node is the depot, and there is no dummy node added. `coordinates` is `n+1` by `2` matrix. 



