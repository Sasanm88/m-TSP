# To use HGSmTSP module without installing as a package
try 
    using .HGSmTSP
catch e 
    include(joinpath(@__DIR__, "../src/HGSmTSP.jl")) 
    using .HGSmTSP
end


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