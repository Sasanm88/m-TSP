using LinearAlgebra 

include(joinpath(@__DIR__, "../src/main.jl")) 


function test_simple()

    n_vehicles = 4
    n_nodes = 100
    coordinates = rand(n_nodes, 2) .* 1000

    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = norm(coordinates[i, :] - coordinates[j, :])
        end
    end

    sol = solve_mTSP(n_vehicles, dist_mtx, coordinates; n_iterations=2000, time_limit=Inf)

end

@time test_simple()