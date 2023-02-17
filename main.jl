include("Create_Sample.jl")
include("MIP_SPLIT.jl")
include("Split.jl")
include("GA.jl")
include("Initial.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Neighborhood.jl")
include("Neighborhood_intra.jl")

function test()
    T = Read_TSPLIB_instance(:berlin52, 1)
    n = size(T)[1]-2
    demands = ones(Int, n)
    K = 5
    W = 35
    h = 0.2
    popsize = (10,50)
    k_tournament = 3
    num_iter = 20000
    Mutation_Chance = 0.1
    num_runs = 10
    avg = 0.0

    P = Chromosome[]
    for i=1:num_runs
        P = Perform_Genetic_Algorithm(T, demands,K, W, h, popsize, k_tournament, num_iter, Mutation_Chance);
        avg += P[1].fitness
    end

    println("Average: ", avg/num_runs)
end

test()