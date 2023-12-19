using Clustering
using Distances
using Random
using StatsBase
using TSPLIB
using TSPSolvers

using LinearAlgebra


include("create_sample.jl")
include("split.jl")
include("genetic_algorithm.jl")
include("initial.jl")
include("mutation.jl")
include("crossover.jl")
include("neighborhood.jl")
include("neighborhood_intra.jl")
include("costs.jl")
include("intersection.jl")
include("enrichment.jl")




"""
    In this function, `dist_mtx` does not have the dummy node added to the end.
"""
function solve_mTSP(
    n_vehicles::Int,
    dist_mtx::Matrix{Float64},
    coordinates::Matrix{Float64};
    n_runs::Int=1,
    n_iterations::Int=100,
    time_limit::Float64=10.0,
    W::Int=1000,
    h::Float64=0.3,
    popsize::Tuple{Int,Int}=(10, 20),
    k_tournament::Int=2,
    mutation_chance::Float64=0.0,
    num_nei::Int=2,
    verbose::Bool=false
)::Tuple{Vector{Vector{Int}},Vector{Float64}}

    if n_vehicles == 1
        dist_mtx_int = round.(Int, dist_mtx)

        tour, tour_len = TSPSolvers.solve_tsp(dist_mtx_int; algorithm="HGS", nbIter=n_iterations, timeLimit=time_limit)
        return [tour], [Float64(tour_len)]
    end


    t0 = time()

    n_nodes = size(dist_mtx)[1]
    @assert size(coordinates, 1) == n_nodes

    depot_coordinates = coordinates[1, :]
    customer_coordinates = coordinates[2:end, :]
    @assert length(depot_coordinates) == 2
    @assert size(customer_coordinates, 1) == n_nodes - 1

    dist_mtx_with_dummy = [
        dist_mtx dist_mtx[:, 1]
        dist_mtx[1, :]' 0.0
    ]

    best = Inf
    worst = 0.0
    crossover_functions = Int[2, 3]
    avg = 0.0

    best_chrm = Chromosome(Int[], 0.0, 0.0, Tour[])
    worst_chrm = Chromosome(Int[], 0.0, 0.0, Tour[])
    all_chrms = Chromosome[]

    for _ in 1:n_runs
        time_limit_for_this_run = time_limit - (time() - t0)
        P, roullet = perform_genetic_algorithm(
            dist_mtx_with_dummy,
            n_vehicles,
            h,
            popsize,
            k_tournament,
            n_iterations,
            time_limit_for_this_run,
            mutation_chance,
            num_nei,
            crossover_functions,
            customer_coordinates,
            depot_coordinates,
            verbose = verbose
        )

        avg += P[1].fitness
        push!(all_chrms, P[1])
        if P[1].fitness < best
            best = P[1].fitness
            best_chrm = P[1]
        end
        if P[1].fitness > worst
            worst = P[1].fitness
            worst_chrm = P[1]
        end
    end

    hgs_routes = [t.sequence for t in best_chrm.tours]
    hgs_route_lengths = [t.cost for t in best_chrm.tours]

    return hgs_routes, hgs_route_lengths
end


function solve_instances_set1(num_customers::Int, num_salesmen::Int, num_instances::Int, verbose::Bool)
    for instance in 1:num_instances
        T, depot, Customers = create_random_sample(num_customers)
        n = size(T)[1] - 2
        h = 0.3
        popsize = (10, 20)
        k_tournament = 2
        num_iter = 2500
        time_limit = Inf
        Mutation_Chance = 0.0
        num_runs = 10
        num_nei = 2
        avg = 0.0
        best = Inf
        worst = 0.0
        crossover_functions = [2]


        P = Chromosome[]
        for run in num_runs
            avg = 0.0
            best = Inf
            worst = 0.0
            t1 = time()
            for i in 1:num_runs
                P, roullet = perform_genetic_algorithm(T, num_salesmen, h, popsize,
                    k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, crossover_functions, Customers, depot, verbose = verbose)
                avg += P[1].fitness
                if P[1].fitness < best
                    best = P[1].fitness
                end
                if P[1].fitness > worst
                    worst = P[1].fitness
                end
            end
            t2 = time()
            println("Results for instance number ", instance, " ,m=", num_salesmen)
            println("Best: ", round(best, digits=1), "  Average: ", round(avg / num_runs, digits=1),
                "  Worst: ", round(worst, digits=1), " , run time= ", round((t2 - t1) / num_runs, digits=1))
        end
    end
end

function solve_instances_set2(instances::Vector{Symbol}, Ms::Vector{Int}, verbose::Bool)
    for instance in instances
        for K in Ms
            T = read_TSPLIB_instance(instance, 1)
            n = size(T)[1] - 2
            tsp = readTSPLIB(instance)
            allNodes = tsp.nodes
            depot = allNodes[1, :]

            Customers = allNodes[2:n+1, :]
            h = 0.3
            popsize = (10, 20)
            k_tournament = 2
            num_iter = 2500
            time_limit = Inf
            Mutation_Chance = 0.0
            num_runs = 10
            num_nei = 2
            avg = 0.0
            best = Inf
            worst = 0.0
            crossover_functions = [2]


            P = Chromosome[]
            for run in num_runs
                avg = 0.0
                best = Inf
                worst = 0.0
                t1 = time()
                for i in 1:num_runs
                    P, roullet = perform_genetic_algorithm(T, K, h, popsize,
                        k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, crossover_functions, Customers, depot, verbose = verbose)
                    avg += P[1].fitness
                    if P[1].fitness < best
                        best = P[1].fitness
                    end
                    if P[1].fitness > worst
                        worst = P[1].fitness
                    end
                end
                t2 = time()
                println("Results for ", instance, " ,m=", K)
                println("Best: ", round(best, digits=1), "  Average: ", round(avg / num_runs, digits=1),
                    "  Worst: ", round(worst, digits=1), " , run time= ", round((t2 - t1) / num_runs, digits=1))
            end
        end
    end
end

function solve_instances_set3_4(dir_name::String, sample_names::Vector{String}, verbose::Bool)

    best_chrm = Chromosome(Int[], 0.0, 0.0, Tour[])
    worst_chrm = Chromosome(Int[], 0.0, 0.0, Tour[])
    all_chrms = Chromosome[]

    row = 0
    for sample_name in sample_names
        row += 1
        println(sample_name)
        m, T, depot_, customers = read_data(dir_name, sample_name)
        Customers = copy(transpose(customers))
        depot = Float64.(depot_)
        n = size(T)[1] - 2
        demands = ones(Int, n)
        W = 1000
        h = 0.1
        popsize = (10, 20)
        k_tournament = 2
        num_iter = 1000000
        time_limit = (n + 1) / 5
        Mutation_Chance = 0.0
        num_runs = 20
        num_nei = 2
        avg = 0.0
        best = Inf
        worst = 0.0
        crossover_functions = Int[2]

        t1 = time()
        for i in 1:num_runs
            P, roullet = perform_genetic_algorithm(T, m, h, popsize,
                k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, crossover_functions, Customers, depot, verbose = verbose)

            avg += P[1].fitness
            push!(all_chrms, P[1])
            if P[1].fitness < best
                best = P[1].fitness
                best_chrm = P[1]
            end
            if P[1].fitness > worst
                worst = P[1].fitness
                worst_chrm = P[1]
            end
        end
        t2 = time()
        println("Results for ", sample_name, " ,m=", m)
        println("Best: ", round(best, digits=2), "  Average: ", round(avg / num_runs, digits=2),
            "  Worst: ", round(worst, digits=2), " , run time= ", round((t2 - t1) / num_runs, digits=0))
    end

    return best_chrm
end





