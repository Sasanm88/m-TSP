This repository is private.

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




## How to compile the package and save the state?

Use [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl).

Create a compiled sysimage by running PackageCompiler as follows:
```julia
using PackageCompiler
create_sysimage(["HGSmTSP"], sysimage_path="sysimage/sys_hgsmtsp.so", precompile_execution_file="sysimage/precompile.jl")
```
This will take for a while. After the compilation is finished, restart Julia by
```
julia --sysimage sysimage/sys_hgsmtsp.so
```
Then test.

You can also do
```
julia --sysimage sysimage/sys_hgsmtsp.so test/test_module.jl
```

## How to the sysimage from Python via pyjulia

First install [`pyjulia`](https://github.com/JuliaPy/pyjulia).

Then use the sysimage as follows:
```python
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
```

NOTE: This is not as fast as expected. Not sure why.
