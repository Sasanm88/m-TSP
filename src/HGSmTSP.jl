module HGSmTSP

    include("main.jl")

    export solve_mTSP




    # for precompilation of the module

    n_nodes, n_vehicles = 10, 2
    coordinates = rand(n_nodes, 2) .* 1000    
    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = norm(coordinates[i, :] - coordinates[j, :])
        end
    end    
    solve_mTSP(n_vehicles, dist_mtx, coordinates)

    
    n_nodes, n_vehicles = 10, 1
    coordinates = rand(n_nodes, 2) .* 1000    
    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = norm(coordinates[i, :] - coordinates[j, :])
        end
    end    
    solve_mTSP(n_vehicles, dist_mtx, coordinates; time_limit=1.0)

end