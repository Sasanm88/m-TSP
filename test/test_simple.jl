using LinearAlgebra 

include(joinpath(@__DIR__, "../src/main.jl")) 

using Random


function test_simple(n_nodes, n_vehicles; n_iterations=1000, time_limit=4.0)
    Random.seed!(1)

    coordinates = rand(n_nodes, 2) .* 10

    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = norm(coordinates[i, :] - coordinates[j, :])
        end
    end

    hgs_routes, hgs_route_lengths = solve_mTSP(n_vehicles, dist_mtx, coordinates; n_iterations=n_iterations, time_limit=time_limit)
    
    return hgs_routes, hgs_route_lengths
end

# @time test_simple(10, 3)

# @time test_simple(10, 1)

@time test_simple(100, 5, time_limit=10000.0)


using BenchmarkTools
@btime test_simple(100, 5, time_limit=10000.0)
# before @view: 1.818 s (25476236 allocations: 2.83 GiB)
# after  @view: 1.333 s (25476236 allocations: 2.28 GiB)

println("Done!")
