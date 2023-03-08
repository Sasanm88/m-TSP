mutable struct Tour
    Sequence::Vector{Int}
    cost::Float64
end

mutable struct Chromosome
    genes::Vector{Int64}
    feasible::Bool
    fitness::Float64      
    power::Float64
    tours::Vector{Tour}
end


function Creat_Random_Cromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i = 1:n_nodes])
    chromosome
end

function Parent_Selection_RWS(Population::Vector{Chromosome}, total_p::Float64, popsize::Int64)  #Random Wheel Selection
    r = -rand() * total_p
    summ = 0
    for i = 1:popsize
        summ -= Population[i].power
        if r > summ
            return Population[i]
        end
    end
end

function Parent_Selection_TS(Population::Vector{Chromosome}, k::Int64, popsize::Int64)  #Tournament Selection
    idx = sample(1:popsize, k, replace=false)
    return Population[idx[argmin(idx)]]
end


function Parent_Selection_RkS(Population::Vector{Chromosome}, ranks::Vector{Int64}, total_r::Int64, popsize::Int64)  #Rank Selection
    r = rand() * total_r
    summ = 0
    for i = 1:popsize
        summ += ranks[i]
        if r < summ
            return Population[i]
        end
    end
end

function Select_parents(Population::Vector{Chromosome}, k_tournament::Int64, popsize::Int64)
    return Parent_Selection_TS(Population, k_tournament, popsize), Parent_Selection_TS(Population, k_tournament, popsize)
end

function Reproduce(TT::Matrix{Float64}, parent1::Chromosome, parent2::Chromosome, n_nodes::Int64)
    #     r = sample(1:length(crsovr_chances), Weights(crsovr_chances),1)
    r = rand(1:2)
#     r = 1
    if r == 10
        return Crossover_OX1(parent1.genes, parent2.genes, n_nodes)  #4730
    elseif r == 20
        return Crossover_OX2(parent1.genes, parent2.genes, n_nodes)  #4618
    elseif r == 30
        return Crossover_POS(parent1.genes, parent2.genes, n_nodes)  #4712
    elseif r == 40
        return Crossover_CX(parent1.genes, parent2.genes, n_nodes)  #4706
    elseif r == 1
        return Crossover_HX(TT, parent1.genes, parent2.genes, n_nodes) #4244
    elseif r==6
        return Crossover_PMX(parent1.genes, parent2.genes, n_nodes)  #4695
    else
        return new_crossover(parent1, parent2, TT, n_nodes)
    end
end

function find_difference(c1::Vector{Int64}, c2::Vector{Int64})  #range between zero and 1, zero when two chromosomes are exactly the same, 1 when all genes are different
    diff1 = 0
    diff2 = 0
    c3 = reverse(c2)
    for i = 1:length(c1)
        if c1[i] != c2[i]
            diff1 += 1
        end
        if c1[i] != c3[i]
            diff2 += 1
        end
    end
    return min(diff1, diff2) / length(c1)
end

function find_neighbors(a::Int, i::Int, m::Int)
    neighbors = Int[]
    for j=max(1, i-m):i-1
        push!(neighbors, j)
    end
    for j=i+1:min(a, i+m)
        push!(neighbors, j)
    end
    return neighbors
end

function Sort_based_on_power(Population::Vector{Chromosome}, num_nei::Int)
    popsize = length(Population)
    diff1 = 0.0
    diff2 = 0.0
    for i = 1:popsize
        neighbors = find_neighbors(popsize, i, num_nei)
        diff = 0.0
        for j in neighbors
            diff1 = find_difference(Population[i].genes, Population[j].genes)
            diff += diff1
        end
        Population[i].power = Population[i].fitness * 0.8^(diff/length(neighbors))
#         if i == 1
#             diff1 = find_difference(Population[1].genes, Population[2].genes)
#             Population[i].power = Population[i].fitness * 0.8^diff1 #(2-diff1)
#         elseif i == popsize
#             Population[i].power = Population[i].fitness * 0.8^diff1 #(2-diff1)
#         else
#             diff2 = find_difference(Population[i].genes, Population[i+1].genes)
#             Population[i].power = Population[i].fitness * 0.8^((diff1 + diff2) / 2) #(2-(diff1+diff2)/2)
#             diff1 = diff2
#         end
    end
    sort!(Population, by=x -> x.power)
end


function Survive(Population::Vector{Chromosome}, mu::Int64, sigma::Int64, count::Int, feasibility::Bool)
    if count < mu+sigma 
        return count
    end

    del_count = 0
    del_idx = Vector{Int64}()
    for i = 1:length(Population)-1
        if del_count == sigma
            break
        end
        if Population[i].feasible == feasibility
            for j = i+1:length(Population)
                if del_count == sigma
                    break
                end
                if Population[j].feasible == feasibility
                    if Population[i].genes == Population[j].genes
                        if !(j in del_idx)
                            push!(del_idx, j)
                            del_count += 1
                        end
                    end
                end
            end
        end
    end
    deleteat!(Population, sort(del_idx))
    del_idx = Vector{Int64}()
    last_index = length(Population)
    index = 0

    while del_count < sigma
        i = last_index - index
        if Population[i].feasible == feasibility
            push!(del_idx, i)
            del_count += 1
        end
        index += 1
    end
    deleteat!(Population, sort(del_idx))
    

    count -= sigma
    return count
end



function Perform_Survival_Plan(Population::Vector{Chromosome}, mu::Int64, sigma::Int64, feas_count::Int, Inf_count::Int)
    feas_count = Survive(Population, mu, sigma, feas_count, true)
    Inf_count = Survive(Population, mu, sigma, Inf_count, false)
    return feas_count, Inf_count
end

function best_route(Population::Vector{Chromosome})
    for tour in Population[1].tours
        for i in tour.Sequence
            print(i, ", ")
        end
        println()
    end
end

function best_objective(Population::Vector{Chromosome})
    for p in Population
        if p.feasible
            return p.fitness
        end
    end
end

function Generate_new_generation(TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, K::Int, W::Int,
        Population::Vector{Chromosome}, popsize::Tuple{Int64,Int64}, k_tournament::Int64, 
        ClosenessT::Matrix{Int64}, Gen_num::Int64, old_best::Float64, improve_count::Int64, Mutation_Chance::Float64,
        tsp_tour::Vector{Int}, roullet::Vector{Int}, num_nei::Int, targetFeasible::Float64, fractionFeasibleLoad::Float64,
        penalty::Float64, max_inf::Int, feas_count::Int, Inf_count::Int)
    t1 = time()

    mu, sigma = popsize
    n_nodes = length(Population[1].genes)

    if fractionFeasibleLoad < targetFeasible - 0.05
        penalty = min(penalty * 1.02, 3.0)
    elseif fractionFeasibleLoad > targetFeasible - 0.05
        penalty = max(penalty * 0.98, 1.2)
    end

    
    if improve_count % 100 == 0  
        Diversify(Population, TT, demands, K, W, mu, tsp_tour, penalty, max_inf)
    end
    Sort_based_on_power(Population, num_nei)
#     println("B")
    psize = length(Population)
    parent1, parent2 = Select_parents(Population, k_tournament, psize)

    child = Reproduce(TT, parent1, parent2, n_nodes)

    offspring, offspring_f = SPLIT_offspring(TT, demands, K, W, child, penalty, max_inf)
    if offspring.feasible
        fractionFeasibleLoad = fractionFeasibleLoad * 0.99 + 0.01
        offspring = Improve_chromosome(offspring, TT, Close_nodes, demands, W, n_nodes, roullet)
        feas_count += 1
    elseif rand() < 0.5
        fractionFeasibleLoad = fractionFeasibleLoad * 0.99   
        feas_count += 1
        offspring = Improve_chromosome(offspring_f, TT, Close_nodes, demands, W, n_nodes, roullet)
    else
        fractionFeasibleLoad = fractionFeasibleLoad * 0.99   #Improve or not?
        Inf_count += 1
    end
    
    push!(Population, offspring)
    sort!(Population, by=x -> x.fitness)
    
    feas_count, Inf_count = Perform_Survival_Plan(Population, mu, sigma, feas_count, Inf_count)

    new_best = best_objective(Population)
    if (old_best - new_best) / new_best > 0.0005
        old_best = new_best
        improve_count = 0
    else
        improve_count += 1
    end
    t2 = time()
    

    if Gen_num % 100 == 0
#         println("Generation ", Gen_num, " the best objective is: ", old_best, "  feasible: ", round(fractionFeasibleLoad, digits = 3), "  penalty: ", round(penalty, digits = 3), "   " , feas_count,"  ", Inf_count)
    end
    Gen_num += 1
    return Gen_num, old_best, Population, improve_count, feas_count, Inf_count, penalty, fractionFeasibleLoad
end

function Perform_Genetic_Algorithm(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int,h::Float64, popsize::Tuple{Int64,Int64},
    k_tournament::Int64, num_iter::Int64, Mutation_Chance::Float64, num_nei::Int, targetFeasible::Float64, max_inf::Int)
    n_nodes = size(TT)[1] - 2
    t1 = time()
    ClosenessT= Find_Closeness(TT, h) 
    mu, sigma = popsize
    penalty = 1.1
    improve_count = 0
    Gen_num = 0
    old_best = 0.0
    roullet = ones(Int, 6) * 100
    tsp_tour = find_tsp_tour1(TT[1:n_nodes+1, 1:n_nodes+1])
    Population, old_best = Generate_initial_population(TT, demands, K, W, mu, tsp_tour, penalty, max_inf) 
    count = 0
    feas_count = mu
    Inf_count = mu
    
    fractionFeasibleLoad = 0.2
    @inbounds while improve_count < num_iter
            Gen_num, old_best, Population, improve_count, feas_count, Inf_count, penalty, fractionFeasibleLoad = Generate_new_generation(TT, ClosenessT, demands, K, W, Population, popsize, k_tournament, ClosenessT, Gen_num, old_best, improve_count, Mutation_Chance, tsp_tour, roullet, num_nei, targetFeasible, fractionFeasibleLoad, penalty, max_inf, feas_count, Inf_count)
        count += 1
    end
    t2 = time()

#     println("The best objective achieved in ", Gen_num, " generations is: ", Population[1].fitness, " and it took ", t2 - t1, " seconds.")
#     println("And the best route is: ")
#     best_route(Population)
    return Population
end


