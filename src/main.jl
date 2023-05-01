include("Create_Sample.jl")
include("Split.jl")
include("GA.jl")
include("Initial.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Neighborhood.jl")
include("Neighborhood_intra.jl")
include("costs.jl")
include("intersection.jl")

function Solve_instances(dir_name::String, sample_names::Vector{String})
    
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
        n = size(T)[1]-2
        demands = ones(Int, n)
        W = 1000
        h = 0.3
        popsize = (10,20)
        k_tournament = 2
        num_iter = 1000000
        time_limit = (n+1)*240/100
        Mutation_Chance = 0.0
        num_runs = 20
        num_nei = 2
        avg = 0.0
        best = Inf
        worst = 0.0
        crossover_functions = [2, 3]

        t1 = time()
        for i=1:num_runs
            P, roullet = Perform_Genetic_Algorithm(T, m, h, popsize, 
                        k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, crossover_functions, Customers, depot);

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
        println("Best: ", round(best, digits = 2), "  Average: ", round(avg/num_runs, digits = 2), 
            "  Worst: ", round(worst, digits = 2), " , run time= ", round((t2-t1)/num_runs, digits=0))
    end
end

function test(instances::Vector{Symbol}, Ms::Vector{Int})
    for instance in instances
        for K in Ms
            T = Read_TSPLIB_instance(instance, 1)
            n = size(T)[1]-2
            tsp = readTSPLIB(instance)
            allNodes = tsp.nodes
            depot = allNodes[1, :]

            Customers = allNodes[2:n+1, :]
            h = 0.3
            popsize = (10,20)
            k_tournament = 2
            num_iter = 2500
            time_limit = Inf
            Mutation_Chance = 0.0
            num_runs = 20
            num_nei = 2
            avg = 0.0
            best = Inf
            worst = 0.0
            crossover_functions = [1, 2, 3]
            
            
            P = Chromosome[]
            for run in num_runs
                avg = 0.0
                best = Inf
                worst = 0.0
                t1 = time() 
                for i=1:num_runs
                    P, roullet = Perform_Genetic_Algorithm(T,K, h, popsize, 
                k_tournament, num_iter, time_limit, Mutation_Chance, num_nei, crossover_functions, Customers, depot);
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

# dir_name = "set1"
# sample_names = ["mtsp150_3", "mtsp150_5", "mtsp150_10", "kroa200_3", "kroa200_5","kroa200_10", "lin318_3", "lin318_5", "lin318_10"]

# Solve_instances(dir_name, sample_names)
    
# instances = [:eil51]
# Ms = [2,3,5,7]
# test(instances, Ms)




