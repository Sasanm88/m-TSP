include("Create_Sample.jl")
include("MIP_SPLIT.jl")
include("Split.jl")
include("GA.jl")
include("Initial.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Neighborhood.jl")
include("Neighborhood_intra.jl")
include("costs.jl")

function test(instances::Vector{Symbol}, Ms::Vector{Int})
    for instance in instances
        for K in Ms
            T = Read_TSPLIB_instance(instance, 1)
            n = size(T)[1]-2
            demands = ones(Int, n)
            W = 150
            h = 0.3
            popsize = (10,50)
            k_tournament = 3
            num_iter = 2000
            Mutation_Chance = 0.1
            num_runs = 10
            
            
            P = Chromosome[]
            Mutation_Chances = [0.0]
            for Mutation_Chance in Mutation_Chances
                avg = 0.0
                best = Inf
                worst = 0.0
                # println("Mutation Chance: ", Mutation_Chance)
                t1 = time() 
                for i=1:num_runs
                    P = Perform_Genetic_Algorithm(T, demands,K, W, h, popsize, k_tournament, num_iter, Mutation_Chance);
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
                println("Best: ", round(best, digits = 1), "  Average: ", round(avg/num_runs, digits = 1), 
                    "  Worst: ", round(worst, digits = 1), " , run time= ", round((t2-t1)/num_runs, digits=1))
            end
        end
    end
end


instances = [:eil51]
Ms = [2,3,5,7]
test(instances, Ms)




