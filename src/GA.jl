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


function Reproduce(TT::Matrix{Float64}, parent1::Chromosome, parent2::Chromosome, n_nodes::Int64, crossover_functions::Vector{Int}) :: Vector{Int}
    r::Int = rand(crossover_functions)
    if r == 1
        return Crossover_HX(TT, parent1.genes, parent2.genes, n_nodes)
    elseif r == 2
        return tour_crossover2(parent1, parent2, TT, n_nodes)
    elseif r == 3
        return tour_crossover3(parent1, parent2, TT, n_nodes)
    elseif r == 4
        return tour_crossover4(parent1, parent2, TT, n_nodes)
    elseif r == 5
        return tour_crossover5(parent1, parent2, TT, n_nodes)
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

function find_difference(c1::Chromosome, c2::Chromosome)  #range between zero and 1, zero when two chromosomes are exactly the same, 1 when all genes are different
    m = length(c1.tours)
    n = length(c1.genes)
    A = zeros(Int, m,m)
    for i=1:m
        for j=1:m
            A[i,j] = length(intersect(Set(c1.tours[i].Sequence), Set(c2.tours[j].Sequence)))
        end
    end
        summ = 0
    while true
        idx = argmax(A)

        if A[idx] == 0
            break
        end
        i = idx[1]
        j = idx[2]
        summ += A[i,j]
        A[i,:] .= 0
        A[:,j] .= 0
    end
    return 1-summ/n
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

function Sort_based_on_power!(Population::Vector{Chromosome}, num_nei::Int)
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


function Perform_Survival_Plan!(Population::Vector{Chromosome}, mu::Int64, sigma::Int64)
    if length(Population) >= mu + sigma
        del_count = 0
        del_idx = Vector{Int64}()
        @inbounds for i = 1:length(Population)-1
            if del_count == length(Population) - mu
                break
            end
            @inbounds for j = i+1:length(Population)
                if del_count == length(Population) - mu
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

function educate_and_add_the_offspring!(offspring::Chromosome, Population::Vector{Chromosome}, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, Customers::Matrix{Float64}, depot::Vector{Float64}, old_best::Float64, roullet::Vector{Int}, n_nodes::Int)
    if n_nodes < 700 && length(offspring.tours) < 20
        if rand() < 0.3
            Solve_all_intersections!(offspring, Customers, depot, TT)
            Enrich_the_chromosome!(offspring, TT, Customers, depot, n_nodes)
        end
        if rand() < 0.3
            Solve_all_intersections!(offspring, Customers, depot, TT)
            Enrich_the_chromosome!(offspring, TT, Customers, depot, n_nodes)
        end
    else
        if rand() < 0.1
            Solve_all_intersections!(offspring, Customers, depot, TT)
        end
    end
    Enrich_the_chromosome!(offspring, TT, Customers, depot, n_nodes)
    Improve_chromosome!(offspring, TT, Close_nodes, n_nodes, roullet, old_best)
        
    push!(Population, offspring)
end

function Generate_new_generation(TT::Matrix{Float64}, Close_nodes::Matrix{Int}, K::Int,
        Population::Vector{Chromosome}, popsize::Tuple{Int64,Int64}, k_tournament::Int64, 
        Gen_num::Int64, old_best::Float64, improve_count::Int64, Mutation_Chance::Float64,
        tsp_tours::Vector{Vector{Int}}, roullet::Vector{Int}, num_nei::Int, crossover_functions::Vector{Int}, Customers::Matrix{Float64}, depot::Vector{Float64}, t0::Float64, time_limit::Float64)
    t1 = time()

    mu, sigma = popsize
    n_nodes = length(Population[1].genes)


    if improve_count % 1000 == 999 
        Diversify!(Population, TT, K, mu, tsp_tours, Customers, depot, improve_count)
    end
    

    Sort_based_on_power!(Population, num_nei)
    psize = length(Population)

    parent1, parent2 = Select_parents(Population, k_tournament, psize)
    child::Vector{Int} = Reproduce(TT, parent1, parent2, n_nodes, crossover_functions)   
    obj, trips = SPLIT(TT, K, child)
    offspring = Chromosome(child, obj, 0.0, trips)
    educate_and_add_the_offspring!(offspring, Population, TT, Close_nodes, Customers, depot, old_best, roullet, n_nodes)

    
    sort!(Population, by=x -> x.fitness)


    Perform_Survival_Plan!(Population, mu, sigma)

    new_best = Population[1].fitness
    if round(old_best, digits=3) > round(new_best, digits=3)
        old_best = new_best
        improve_count = 0
    else
        improve_count += 1
    end
    t2 = time()
    

    if Gen_num % 10 == 0
        println("Generation ", Gen_num, " the best objective is: ", old_best, "   time left: $(round(t0+time_limit -time())) seconds")
    end
    Gen_num += 1
    return Gen_num, old_best, Population, improve_count
end


function Perform_Genetic_Algorithm(
    TT::Matrix{Float64}, K::Int, h::Float64, popsize::Tuple{Int64,Int64},
    k_tournament::Int64, num_iter::Int64, time_limit::Float64, Mutation_Chance::Float64, num_nei::Int, crossover_functions::Vector{Int},
    Customers::Matrix{Float64}, depot::Vector{Float64}
) 
    Random.seed!(Int(round(time())))

    n_nodes = size(TT)[1] - 2
    t1 = time()
    ClosenessT= Find_Closeness(TT, h) 
    mu, sigma = popsize
    improve_count = 0
    Gen_num = 0
    old_best = 0.0
    roullet = ones(Int, 6) * 100
    
    tsp_tours = find_tsp_tour2(TT[1:n_nodes+1, 1:n_nodes+1])
    
    if n_nodes < 1400
        tsp_tour, _ = find_tsp_tour1(TT[1:n_nodes+1, 1:n_nodes+1])
        push!(tsp_tours, tsp_tour)
    end
    
    Population, old_best = Generate_initial_population(TT, K, mu, tsp_tours, Customers, depot) 

    count = 0

    while improve_count < num_iter
        if time() - t1 >= time_limit
            break
        end

        Gen_num, old_best, Population, improve_count = Generate_new_generation(TT, ClosenessT, K,
        Population, popsize, k_tournament, Gen_num, old_best, improve_count, Mutation_Chance, tsp_tours, roullet, num_nei, crossover_functions, Customers, depot, t1, time_limit)

        count += 1
    end
    t2 = time()

    println("The best objective achieved in ", Gen_num, " generations is: ", Population[1].fitness, " and it took ", t2 - t1, " seconds.")
#     println("And the best route is: ")
#     best_route(Population)
    return Population, roullet
end


