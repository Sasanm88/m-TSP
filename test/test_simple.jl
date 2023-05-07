using LinearAlgebra 

include(joinpath(@__DIR__, "../src/main.jl")) 


function test_simple(n_nodes, n_vehicles; n_iterations=10000, time_limit=4.0)

    coordinates = rand(n_nodes, 2) .* 10

    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = norm(coordinates[i, :] - coordinates[j, :])
        end
    end

    hgs_routes, hgs_route_lengths = solve_mTSP(n_vehicles, dist_mtx, coordinates; n_iterations=n_iterations, time_limit=time_limit)

end

@time test_simple(10, 3)

@time test_simple(10, 1)

@time test_simple(100, 5)