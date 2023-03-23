mutable struct Tour
    Sequence::Vector{Int}
    cost::Float64
end

mutable struct Chromosome
    genes::Vector{Int64}
    fitness::Float64
#     feasible::Bool        #We will consider the infeasible chromosomes later
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

function Reproduce(TT::Matrix{Float64}, parent1::Chromosome, parent2::Chromosome, n_nodes::Int64, crossover_functions::Vector{Int})
    #     r = sample(1:length(crsovr_chances), Weights(crsovr_chances),1)
    rs = crossover_functions
    r = rs[rand(1:length(rs))]
    
#     r = 7
    if r == 1
        return Crossover_OX1(parent1.genes, parent2.genes, n_nodes), r  #4730
    elseif r == 2
        return Crossover_OX2(parent1.genes, parent2.genes, n_nodes), r  #4618
    elseif r == 3
        return Crossover_POS(parent1.genes, parent2.genes, n_nodes), r  #4712
    elseif r == 4
        return Crossover_CX(parent1.genes, parent2.genes, n_nodes), r  #4706
    elseif r == 5
        return Crossover_HX(TT, parent1.genes, parent2.genes, n_nodes), r #4244
    elseif r == 6
        return Crossover_PMX(parent1.genes, parent2.genes, n_nodes), r  #4695
    elseif r == 7
        return new_crossover(parent1, parent2, TT, n_nodes), r
    elseif r == 8
        return tour_crossover(parent1, parent2, TT, n_nodes), r
    elseif r == 9
        return tour_crossover2(parent1, parent2, TT, n_nodes), r
    elseif r == 10
        return tour_crossover3(parent1, parent2, TT, n_nodes), r
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

    end
    sort!(Population, by=x -> x.power)
end


function Perform_Survival_Plan(Population::Vector{Chromosome}, mu::Int64, sigma::Int64)
    if length(Population) >= mu + sigma
        del_count = 0
        del_idx = Vector{Int64}()
        @inbounds for i = 1:length(Population)-1
            if del_count == sigma
                break
            end
            @inbounds for j = i+1:length(Population)
                if del_count == sigma
                    break
                end
                if Population[i].genes == Population[j].genes
                    if !(j in del_idx)
                        push!(del_idx, j)
                        del_count += 1
                    end
                end
            end
        end
        deleteat!(Population, sort(del_idx))
        del_idx = Vector{Int64}()
        last_index = length(Population)
        index = 0

        @inbounds while del_count < sigma
            i = last_index - index
            push!(del_idx, i)
            del_count += 1
            index += 1
        end
        deleteat!(Population, sort(del_idx))
    end
end


function best_route(Population::Vector{Chromosome})
    for tour in Population[1].tours
        for i in tour.Sequence
            print(i, ", ")
        end
        println()
    end
end

function Generate_new_generation(TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, K::Int, W::Int,
        Population::Vector{Chromosome}, popsize::Tuple{Int64,Int64}, k_tournament::Int64, 
        ClosenessT::Matrix{Int64}, Gen_num::Int64, old_best::Float64, improve_count::Int64, Mutation_Chance::Float64,
        tsp_tour::Vector{Int}, roullet::Vector{Int}, num_nei::Int, crossover_functions::Vector{Int})
    t1 = time()

    mu, sigma = popsize
    n_nodes = length(Population[1].genes)

    if improve_count % 200 == 19  
        Diversify(Population, TT, demands, K, W, mu, tsp_tour)
    end
    
    Sort_based_on_power(Population, num_nei)
    psize = length(Population)
    parent1, parent2 = Select_parents(Population, k_tournament, psize)

    child, crss = Reproduce(TT, parent1, parent2, n_nodes, crossover_functions)
    obj, trips = SPLIT(TT, demands, K, W, child)
    offspring = Chromosome(child, obj, 0.0, trips)

#     Mutate(child, Mutation_Chance)
    
    
#     if round(obj,digits=4) < round(old_best, digits=4)
# #         println("Improvement in generation ", Gen_num, " by crossover: " , crss,"  ", round(old_best, digits=4) ," to ", round(obj,digits=4))
#         old_best = obj
#         counter2[crss] += 1
#     end
    
    if rand() < Mutation_Chance
        if rand() < 1.0
#             prob_mutation(offspring, T, n_nodes, 0.5)
            chunk_length = rand(2:max(5, minimum([length(tour.Sequence) for tour in offspring.tours])))
            offspring = chunk_mutation(offspring, T, chunk_length, n_nodes)
        else
            pm = 0.5*(1 - improve_count/3000)
            offspring = new_mutation(offspring, T, pm)
        end
    end
    offspring, imprv = Improve_chromosome(offspring, TT, Close_nodes, demands, W, n_nodes, roullet)
    
#     if round(offspring.fitness, digits=4) < round(old_best, digits=4)
# #         println("Improvement in generation ", Gen_num, " by neighborhoods: ", imprv ,"  ", round(old_best, digits=4) ," to ", round(offspring.fitness, digits=4))
#         old_best = offspring.fitness
#         counter1[imprv] += 1
#     end
    
    push!(Population, offspring)
    sort!(Population, by=x -> x.fitness)

    Perform_Survival_Plan(Population, mu, sigma)
    if improve_count % 200 == 0
        Improve_Population(Population, TT, Close_nodes, demands, W, n_nodes)
    end
        
    new_best = Population[1].fitness
    if (old_best - new_best) / new_best > 0.0005
        old_best = new_best
        improve_count = 0
    else
        improve_count += 1
    end
    t2 = time()
    

#     if Gen_num % 100 == 0
#         println("Generation ", Gen_num, " the best objective is: ", old_best)
#     end
    Gen_num += 1
    return Gen_num, old_best, Population, improve_count
end

function Perform_Genetic_Algorithm(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int,h::Float64, popsize::Tuple{Int64,Int64},
    k_tournament::Int64, num_iter::Int64, Mutation_Chance::Float64, num_nei::Int, crossover_functions::Vector{Int})
    n_nodes = size(TT)[1] - 2
    t1 = time()
    ClosenessT= Find_Closeness(TT, h) 
    mu, sigma = popsize
    improve_count = 0
    Gen_num = 0
    old_best = 0.0
    roullet = ones(Int, 6) * 100
    tsp_tour, _ = find_tsp_tour1(TT[1:n_nodes+1, 1:n_nodes+1])
    Population, old_best = Generate_initial_population(TT, demands, K, W, mu, tsp_tour) 
    count = 0

    @inbounds while improve_count < num_iter
            Gen_num, old_best, Population, improve_count = Generate_new_generation(TT, ClosenessT, demands, K, W,
        Population, popsize, k_tournament, ClosenessT, Gen_num, old_best, improve_count, Mutation_Chance, tsp_tour, roullet, num_nei, crossover_functions)
        count += 1
    end
    t2 = time()

#     println("The best objective achieved in ", Gen_num, " generations is: ", Population[1].fitness, " and it took ", t2 - t1, " seconds.")
#     println("And the best route is: ")
#     best_route(Population)
#     println(roullet)
    return Population, roullet
end


