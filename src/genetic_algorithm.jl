mutable struct Tour
    sequence::Vector{Int}
    cost::Float64
end

mutable struct Chromosome
    genes::Vector{Int64}
    fitness::Float64
    #     feasible::Bool        #We will consider the infeasible chromosomes later
    power::Float64
    tours::Vector{Tour}
end


function parent_selection_RWS(population::Vector{Chromosome}, total_p::Float64, popsize::Int64)  #Random Wheel Selection
    r = -rand() * total_p
    summ = 0
    for i in 1:popsize
        summ -= population[i].power
        if r > summ
            return population[i]
        end
    end
end

function parent_selection_TS(population::Vector{Chromosome}, k::Int64, popsize::Int64)  #Tournament Selection
    idx = sample(1:popsize, k, replace=false)
    return population[idx[argmin(idx)]]
end


function parent_selection_RkS(population::Vector{Chromosome}, ranks::Vector{Int64}, total_r::Int64, popsize::Int64)  #Rank Selection
    r = rand() * total_r
    summ = 0
    for i in 1:popsize
        summ += ranks[i]
        if r < summ
            return population[i]
        end
    end
end

function select_parents(population::Vector{Chromosome}, k_tournament::Int64, popsize::Int64)
    return parent_selection_TS(population, k_tournament, popsize), parent_selection_TS(population, k_tournament, popsize)
end


function reproduce(TT::Matrix{Float64}, parent1::Chromosome, parent2::Chromosome, n_nodes::Int64, crossover_functions::Vector{Int})::Vector{Int}
    r::Int = rand(crossover_functions)
    if r == 1
        return crossover_HX(TT, parent1.genes, parent2.genes, n_nodes)
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
    for i in 1:length(c1)
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
    A = zeros(Int, m, m)
    for i in 1:m
        for j in 1:m
            A[i, j] = length(intersect(Set(c1.tours[i].sequence), Set(c2.tours[j].sequence)))
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
        summ += A[i, j]
        A[i, :] .= 0
        A[:, j] .= 0
    end
    return 1 - summ / n
end

function find_neighbors(a::Int, i::Int, m::Int)
    neighbors = Int[]
    for j = max(1, i - m):i-1
        push!(neighbors, j)
    end
    for j = i+1:min(a, i + m)
        push!(neighbors, j)
    end
    return neighbors
end

function sort_based_on_power!(population::Vector{Chromosome}, num_nei::Int)
    popsize = length(population)
    diff1 = 0.0
    diff2 = 0.0
    for i in 1:popsize
        neighbors = find_neighbors(popsize, i, num_nei)
        diff = 0.0
        for j in neighbors
            diff1 = find_difference(population[i].genes, population[j].genes)
            diff += diff1
        end
        population[i].power = population[i].fitness * 0.8^(diff / length(neighbors))

    end
    sort!(population, by=x -> x.power)
end


function perform_survival_plan!(population::Vector{Chromosome}, mu::Int64, sigma::Int64)
    if length(population) >= mu + sigma
        del_count = 0
        del_idx = Vector{Int64}()
        @inbounds for i in 1:length(population)-1
            if del_count == length(population) - mu
                break
            end
            @inbounds for j = i+1:length(population)
                if del_count == length(population) - mu
                    break
                end
                if population[i].genes == population[j].genes
                    if !(j in del_idx)
                        push!(del_idx, j)
                        del_count += 1
                    end
                end
            end
        end
        deleteat!(population, sort(del_idx))
        del_idx = Vector{Int64}()
        last_index = length(population)
        index = 0

        @inbounds while del_count < sigma
            i = last_index - index
            push!(del_idx, i)
            del_count += 1
            index += 1
        end
        deleteat!(population, sort(del_idx))
    end
end


function best_route(population::Vector{Chromosome})
    for tour in population[1].tours
        for i in tour.sequence
            print(i, ", ")
        end
        println()
    end
end

function educate_and_add_the_offspring!(offspring::Chromosome, population::Vector{Chromosome}, TT::Matrix{Float64}, Close_nodes::Matrix{Bool}, Customers::Matrix{Float64}, depot::Vector{Float64}, old_best::Float64, roullet::Vector{Int}, 
        n_nodes::Int, improve_count::Int)
    if n_nodes > 700 && length(offspring.tours) < 10
        if rand() < 0.1 && improve_count > 100
            solve_all_intersections!(offspring, Customers, depot, TT)
        end
    else
        if rand() < 0.3 
            solve_all_intersections!(offspring, Customers, depot, TT)
            enrich_the_chromosome2!(offspring, TT, Customers, depot, n_nodes)
        end
        if rand() < 0.3 
            solve_all_intersections!(offspring, Customers, depot, TT)
            enrich_the_chromosome2!(offspring, TT, Customers, depot, n_nodes)
        end
    end
    enrich_the_chromosome2!(offspring, TT, Customers, depot, n_nodes)
    if improve_count > 100
        Improve_chromosome!(offspring, TT, Close_nodes, n_nodes, roullet, old_best, 1000)
    else
        Improve_chromosome!(offspring, TT, Close_nodes, n_nodes, roullet, old_best, 100)
    end
    push!(population, offspring)
end

function generate_new_generation(TT::Matrix{Float64}, close_nodes::Matrix{Bool}, K::Int,
    population::Vector{Chromosome}, popsize::Tuple{Int64,Int64}, k_tournament::Int64,
    gen_num::Int64, old_best::Float64, improve_count::Int64, mutation_chance::Float64,
    tsp_tours::Vector{Vector{Int}}, roullet::Vector{Int}, num_nei::Int, crossover_functions::Vector{Int},
    customers::Matrix{Float64}, depot::Vector{Float64}, t0::Float64, time_limit::Float64;
    verbose::Bool=false)
    t1 = time()

    mu, sigma = popsize
    n_nodes = length(population[1].genes)


    if improve_count % 1000 == 999
        diversify!(population, TT, K, mu, tsp_tours, customers, depot, improve_count)
    end

    sort_based_on_power!(population, num_nei)
    psize = length(population)
    
    if rand() < mutation_chance
        offspring = mutate(rand(population[1:5]), TT, n_nodes)
    else
        parent1, parent2 = select_parents(population, k_tournament, psize)
        child::Vector{Int} = reproduce(TT, parent1, parent2, n_nodes, crossover_functions)
        obj, trips = SPLIT(TT, K, child)
        offspring = Chromosome(child, obj, 0.0, trips)
    end
    
    educate_and_add_the_offspring!(offspring, population, TT, close_nodes, customers, depot, old_best, roullet, n_nodes, improve_count)

    sort!(population, by=x -> x.fitness)

    perform_survival_plan!(population, mu, sigma)

    new_best = population[1].fitness
    if round(old_best, digits=3) > round(new_best, digits=3)
        old_best = new_best
        improve_count = 0
    else
        improve_count += 1
    end
    t2 = time()


    if verbose
        if gen_num % 1000 == 0
            println("Generation ", gen_num, " the best objective is: ", old_best, "   time left: $(round(t0+time_limit -time())) seconds")
        end
    end
    gen_num += 1
    return gen_num, old_best, population, improve_count
end


function perform_genetic_algorithm(
    TT::Matrix{Float64}, K::Int, h::Float64, popsize::Tuple{Int64,Int64},
    k_tournament::Int64, num_iter::Int64, time_limit::Float64, mutation_chance::Float64, num_nei::Int, crossover_functions::Vector{Int},
    customers::Matrix{Float64}, depot::Vector{Float64};
    verbose=false
)
    # Random.seed!(Int(round(time())))

    n_nodes = size(TT)[1] - 2
    t1 = time()
    ClosenessT = find_closeness(TT, h)
    mu, sigma = popsize
    improve_count = 0
    Gen_num = 0
    old_best = 0.0
    roullet = ones(Int, 4) * 100

    tsp_tours = find_tsp_tour2(TT[1:n_nodes+1, 1:n_nodes+1])


    if n_nodes < 1200 
        tsp_tour, _ = find_tsp_tour1(TT[1:n_nodes+1, 1:n_nodes+1])
        push!(tsp_tours, tsp_tour)
    end


    Population, old_best = generate_initial_population(TT, K, mu, tsp_tours, customers, depot)

    count = 0

    if verbose
        println("The initialization took ", time() - t1, " seconds.")
    end

    while improve_count < num_iter
        if time() - t1 >= time_limit
            break
        end

        Gen_num, old_best, Population, improve_count = generate_new_generation(TT, ClosenessT, K,
            Population, popsize, k_tournament, Gen_num, old_best, improve_count, mutation_chance, tsp_tours, roullet, num_nei, crossover_functions, customers, depot, t1, time_limit, verbose=verbose)

        count += 1
    end
    t2 = time()

    if verbose 
        println("The best objective achieved in ", Gen_num, " generations is: ", Population[1].fitness, " and it took ", t2 - t1, " seconds.")
        println("And the best route is: ")
        best_route(Population)
    end

    return Population, roullet
end


