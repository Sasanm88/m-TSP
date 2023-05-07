module HGSmTSP

    include("main.jl")

    export solve_mTSP




    # for precompilation of the module
    n_vehicles = 3
    n_nodes = 10
    coordinates = rand(n_nodes, 2) .* 1000

    dist_mtx = zeros(n_nodes, n_nodes)
    for i=1:n_nodes
        for j=1:n_nodes
            dist_mtx[i, j] = sum(abs.(coordinates[i, :] - coordinates[j, :]))
        end
    end

    t0 = time()
    solve_mTSP(n_vehicles, dist_mtx, coordinates)
    println(time() - t0)

end