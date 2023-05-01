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
include("Draw.jl")
include("Escape.jl")
using XLSX

function Write_to_excel(exfile, sheetnumber, instance_count, instance_name,m, best, avg, worst, run_time)
    instance_count += 1
    XLSX.openxlsx(exfile, mode="rw") do xf
        sheet = xf[sheetnumber]
        sheet["A"*string(instance_count)] = string(instance_name)
        sheet["B"*string(instance_count)] = m
        sheet["C"*string(instance_count)] = best
        sheet["D"*string(instance_count)] = avg
        sheet["E"*string(instance_count)] = worst
        sheet["F"*string(instance_count)] = run_time
    end
end

function test_parameters(mu::Int, sigma::Int, h::Float64, k_tournament::Int, num_nei::Int, crossover_functions::Vector{Int}, sheetnumber::Int)
    instances = [:eil51, :berlin52, :eil76, :rat99]
    LKH3 = [[222.7, 159.6, 124.0, 112.1],[4110.2, 3184.2, 2440.9, 2440.9], [280.9, 196.7, 143.4, 128.2],[690.8, 523.3, 467.0,442.5]]
    Ms = [2, 3, 5, 7]

    count = 0
    best_ = 0.0
    Avg_ = 0.0
    worst_ = 0.0
    time_ = 0.0
    P = Chromosome[]
    instance_count = 1
    for (i,instance) in enumerate(instances)
        for (j,K) in enumerate(Ms)
            
            count += 1
            T = Read_TSPLIB_instance(instance, 1)
            n = size(T)[1]-2
            demands = ones(Int, n)
            W = 150
            popsize = (mu,sigma)
            num_iter = 2000
            Mutation_Chance = 0.0
            num_runs = 10
            avg = 0.0
            best = Inf
            worst = 0.0
            t1 = time() 
            P = Chromosome[]
            for i=1:num_runs
                P, roullet = Perform_Genetic_Algorithm(T, demands,K, W, h, popsize, 
                    k_tournament, num_iter, Mutation_Chance, num_nei, crossover_functions);
    #             roullet_ = roullet_ + roullet
                avg += P[1].fitness
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
            println("Results for ", instance, " ,m=", K)
            println("Best: ", round(best, digits = 1), "  Average: ", round(avg/num_runs, digits = 1), 
                "  Worst: ", round(worst, digits = 1), " , run time= ", round((t2-t1)/num_runs, digits=0))
            best_ += 100*(best-LKH3[i][j])/LKH3[i][j]
            worst_ += 100*(worst-LKH3[i][j])/LKH3[i][j]
            Avg_ += 100*(avg/num_runs-LKH3[i][j])/LKH3[i][j]
            time_ += t2-t1
            Write_to_excel("crossover2.xlsx", sheetnumber, instance_count, instance, K, best, avg/num_runs, worst, (t2-t1)/num_runs)
            instance_count += 1
        end
    end
    Write_to_excel("crossover2.xlsx", sheetnumber, instance_count, "summary", "", best_/count, Avg_/count, worst_/count, time_/count)
    println("Results:  Best: ", round(best_/count, digits=1), "  Average: ", round(Avg_/count, digits=1), "  Worst: ", round(worst_/count, digits=1), "   run time: ", round(time_/count, digits=1))
end
function test()
    h = 0.3
    k_tournament = 2
    crossover_functions = [5,7,8] #[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]]
    sheetnumber = 8
    mu = 20
    sigma = 40
    num_nei = 2
    
    test_parameters(mu, sigma, h, k_tournament, num_nei, crossover_functions, sheetnumber)

end

test()