using LKH

using Random, Distances, Clustering


# Function to move the farthest median closer to the depot
function move_farthest_median(medians, depot)
    distances = [euclidean(median, depot) for median in medians]
    farthest_median = argmax(distances)
    medians[farthest_median] += 0.1 * (depot - medians[farthest_median])
    return medians
end


function k_median(Customers::Matrix{Float64}, depot::Vector{Float64}, k::Int)
    data = Customers
    n = size(data)[2]
    # Set the number of clusters you want to create and the depot location
    # depot = [0.5, 0.5]

    # Choose k initial cluster centers (medians)
    medians = [data[:,i] for i in sample(1:n, k, replace = false)]

    # Set a threshold for the difference between the maximum and minimum distances
    # from the cluster medians to the depot
    threshold = 0.1
    assignments_ = Int[]
    # Run the modified k-Median clustering algorithm
    for i in 1:100
        # Assign each node to the nearest median
        assignments_ = [argmin([euclidean(data[:,p], m) for m in medians]) for p in 1:n]

        # Compute the median of each cluster
        for j in 1:k
            cluster_points = data[:,findall(x->x==j, assignments_)]
            medians[j] = median(cluster_points, dims=2)[:,1]
        end

        # Calculate the sum of distances from each cluster median to the depot
        distances = [euclidean(median, depot) for median in medians]
        max_distance = maximum(distances)
        min_distance = minimum(distances)

        # If the difference between the maximum and minimum distances is greater than
        # the threshold, move the farthest median closer to the depot
        if max_distance - min_distance > threshold
            medians = move_farthest_median(medians, depot)
        end
    end
    assignments_
end

function initial_kmedian_solution(T::Matrix{Float64}, Customers::Matrix{Float64}, depot::Vector{Float64}, K::Int)
#     assignments_ = k_median(Customers, depot, K)
    k1, k2 = size(Customers)
    if k2 == 2
        Customers_ = transpose(Customers)
    else 
        Customers_ = Customers
    end
    if rand() < 1
        result = kmeans(Customers_, K)
        assignments_ = copy(result.assignments)
    else
        assignments_ = k_median(Customers, depot, K)
    end
    tours = Tour[]
    genes = Int[]
    obj = 0.0
    for i=1:K
        t1 = findall(x->x==i, assignments_)
        if length(t1) == 1
            obj1 = T[1, t1[1]+1] + T[t1[1]+1, 1]
            push!(tours, Tour(t1, obj1))
            push!(genes, t1[1])
            if obj1 > obj
                obj = obj1
            end
        else
            t2 = copy(t1)
            pushfirst!(t2, 0)
            t2 = t2.+1
            TT = T[t2, t2]
            tt1 , obj1 = find_tsp_tour1(TT)
            push!(tours, Tour(t1[tt1], obj1))
            genes = vcat(genes, t1[tt1])
            if obj1 > obj
                obj = obj1
            end
        end
    end
    
    chrm = Chromosome(genes, obj, 0.0, tours)
#     if rand() < 0.5
#         Solve_all_intersections(chrm, Customers, depot, T)
#     end
    return chrm
end

function initial_random_solution(T::Matrix{Float64}, K::Int, n_nodes::Int)
#     assignments_ = k_median(Customers, depot, K)
    
    assignments_ = rand(1:K, n_nodes-K)
    for i=1:K
        insert!(assignments_, rand(1:n_nodes-K+i), i)
    end
    tours = Tour[]
    genes = Int[]
    obj = 0.0
    for i=1:K
        t1 = findall(x->x==i, assignments_)
        if length(t1) == 1
            obj1 = T[1, t1[1]+1] + T[t1[1]+1, 1]
            push!(tours, Tour(t1, obj1))
            push!(genes, t1[1])
            if obj1 > obj
                obj = obj1
            end
        else
            t2 = copy(t1)
            pushfirst!(t2, 0)
            t2 = t2.+1
            TT = T[t2, t2]
            tt1 , obj1 = find_tsp_tour1(TT)
            push!(tours, Tour(t1[tt1], obj1))
            genes = vcat(genes, t1[tt1])
            if obj1 > obj
                obj = obj1
            end
        end
    end
    chrm = Chromosome(genes, obj, 0.0, tours)
    return chrm
end

function find_tsp_tour1(Ct::Matrix{Float64})
    scale_factor = 1000
    dist_mtx = round.(Int, Ct .* scale_factor)

    # tsp_tour, tsp_tour_len = Concorde.solve_tsp(dist_mtx)
    tsp_tour, tour_length = LKH.solve_tsp(dist_mtx)

    @assert tsp_tour[1] == 1

    return tsp_tour[2:length(tsp_tour)].-1, tour_length/scale_factor
end

function Change_initial(c::Vector{Int64}, n_nodes::Int64)
    cc = copy(c)
    r = rand()
    if r < 0.5
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
        if idx1 > idx2
            temp = idx1
            idx1 = idx2
            idx2 = temp
        end
        cc[idx1:idx2] = reverse(cc[idx1:idx2])
    else
        idx = sample(1:n_nodes, rand(2:Int(floor(n_nodes/2))), replace = false) 
        cc[idx] = shuffle(cc[idx])
    end
    return cc
end

function Creat_Random_Cromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i = 1:n_nodes])
    chromosome
end
function Creat_Random_Cromosome2(T::Matrix{Float64}, n_nodes::Int64, m::Int)
    Nodes = shuffle([i for i=1:n_nodes]);
    tours = Tour[]
    for i=1:m
        push!(tours, Tour(Int[], 0.0))
    end

    while length(Nodes) > 0
        city = pop!(Nodes)
        put_city_in_tour(tours, city, T, n_nodes)
    end

    S = Int[]
    obj = 0.0
    for tour in tours
        S = vcat(S, tour.Sequence)
    end
    return S
end

function Creat_Random_Cromosome3(T::Matrix{Float64}, n_nodes::Int64)
    Nodes = [i for i=1:n_nodes];
    tours = Tour[]
    a = copy(T[1, 2:n_nodes+1])
    b = sortperm(a)

    for i=1:m
        push!(tours, Tour([b[i]], 0.0))
    end

    tour_indices = [i for i=1:m]
    while length(Nodes) > 0
        r = tour_indices[rand(1:length(tour_indices))]
        tour = tours[r]
        last_city = tour.Sequence[length(tour.Sequence)]
        a = copy(T[last_city+1, Nodes.+1])
        new_city = Nodes[argmin(a)]
        push!(tour.Sequence, new_city)
        deleteat!(tour_indices, findfirst(x->x==r, tour_indices))
        deleteat!(Nodes, findfirst(x->x==new_city, Nodes))
        if length(tour_indices) == 0
            tour_indices = [i for i=1:m]
        end
    end

    S = Int[]
    for tour in tours
        S = vcat(S, tour.Sequence)
    end
    return S
end

function Generate_initial_population(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int, tsp_tour::Vector{Int}, Customers::Matrix{Float64}, depot::Vector{Float64})
    Population = Chromosome[]
    n_nodes = length(demands)
    obj, trips = SPLIT(TT, demands, K, W, tsp_tour)
    push!(Population, Chromosome(tsp_tour, obj, 0.0, trips))
    S = Int[]
    for i=1:mu-1
#         print(i, "  ")
        if rand() < 1
#             if rand() < 0.5
            S = Change_initial(tsp_tour, n_nodes)
            obj, trips = SPLIT(TT, demands, K, W, S)
            push!(Population, Chromosome(S, obj, 0.0, trips))
        else
            chrm = initial_kmedian_solution(TT, Customers, depot, K)
            push!(Population, chrm)
        end
        
    end
#     sort!(Population, by=x -> x.fitness)
    return Population, Population[1].fitness
end

function Diversify(Population::Vector{Chromosome}, TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int, tsp_tour::Vector{Int}, Customers::Matrix{Float64}, depot::Vector{Float64}, num::Int)
    n_nodes = length(demands)
    n_best = Int(round(0.15 * mu)) 
#     p = min(0.5, 1 - num/5000)
#     for i=n_best+1:length(Population)
#         if rand() < p
#             S = Change_initial(tsp_tour, n_nodes)
#         else
#             if rand() < 0.5
#                 S = Creat_Random_Cromosome2(TT, n_nodes)
#             else
#                 S = Creat_Random_Cromosome(n_nodes)
#             end
#         end
#         obj, trips = SPLIT(TT, demands, K, W, S)
#         Population[i] = Chromosome(S, obj, 0.0, trips)
#     end
    for i=n_best+1:length(Population)
        S = Creat_Random_Cromosome2(TT, n_nodes, K)
        obj, trips = SPLIT(TT, demands, K, W, S)
        chrm = Chromosome(S, obj, 0.0, trips)
        if rand() < 0.4
            Solve_all_intersections(chrm, Customers, depot, TT)
        end
        if rand() < 0.4
            chrm = Enrich_the_chromosome(chrm, TT, Customers, n_nodes)
        end
        Population[i] = chrm
    end
    sort!(Population, by=x -> x.fitness)
end

function Diversify_(Population::Vector{Chromosome}, TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int, tsp_tour::Vector{Int}, Customers::Matrix{Float64}, depot::Vector{Float64})
    n_nodes = length(demands)
    n_best = Int(round(0.15 * mu)) 
#     println(length(Population))
    for i=n_best+1:length(Population)
#         if rand() < 0.5
# #             S = Change_initial(tsp_tour, n_nodes)
# #             obj, trips = SPLIT(TT, demands, K, W, S)
# #             Population[i] = Chromosome(S, obj, 0.0, trips)
# #         elseif rand() < 0.6
# #             S = Creat_Random_Cromosome(n_nodes)
# #             obj, trips = SPLIT(TT, demands, K, W, S)
# #             Population[i] = Chromosome(S, obj, 0.0, trips)
#             chrm = initial_random_solution(TT, K, n_nodes)
#             Population[i] = chrm
#         else
        chrm = initial_kmedian_solution(TT, Customers, depot, K)
#         chrm = initial_random_solution(TT, K, n_nodes)
        Population[i] = chrm
#         end
    end
    sort!(Population, by=x -> x.fitness)
end

function find_tour_neighbors(tours::Vector{Tour}, Customers::Matrix{Float64}, m::Int)
    means = [mean(Customers[t1, :], dims=1)[1,:] for t1 in [tours[i].Sequence for i=1:m]]
    distances = zeros(m, m)
    for i = 1:m-1
        for j = i+1:m
            distances[i,j] = euclidean(means[i], means[j])
            distances[j,i] = distances[i,j]
        end
    end
    return [sortperm(distances[i,:])[2:min(m,3)] for i=1:m]
end


function Enrich_the_chromosome(Chrm::Chromosome, T::Matrix{Float64}, Customers::Matrix{Float64}, n_nodes::Int)   #Shift(0,1)
    m = length(Chrm.tours)
    max_tour_index = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    max_tour_length = Chrm.fitness
    tour_neighbors = find_tour_neighbors(Chrm.tours, Customers, m)
    improved = true
    count = 0
    while improved && count < 100
        count += 1
        improved = false
        for r1 = 1:m
            for r2 in tour_neighbors[r1]
                if r2 != max_tour_index
                    tour1 = Chrm.tours[r1].Sequence
                    tour2 = Chrm.tours[r2].Sequence
                    
                    k1 = 0
                    became_max = false
                    while k1 < length(tour1) && !became_max
                        cost1 = Chrm.tours[r1].cost
                        cost2 = Chrm.tours[r2].cost
                        k1 += 1
                        city1 = tour1[k1]
                        for k2 = 1:length(tour2)+1
                            new_cost2 = Calculate_new_cost_add_one(tour2, cost2, city1, k2, T, n_nodes)
                            new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, T, n_nodes)
                            do_it = false
                            if r1 == max_tour_index
                                if new_cost2 < max_tour_length
                                    do_it = true
                                end
                            else
                                if (new_cost2 - cost2) < (cost1 - new_cost1) && new_cost2< max_tour_length
                                    do_it = true
                                end
                            end
                            if do_it
                                
                                Chrm.tours[r1].cost = new_cost1
                                Chrm.tours[r2].cost = new_cost2
                                insert!(tour2, k2, city1)  
                                deleteat!(tour1, k1)
                                if new_cost2 > max_tour_length
                                    max_tour_length = new_cost2
                                    max_tour_index = r2 
                                    became_max = true
                                    
                                end
                                k1 -= 1
                                improved = true
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    for tour in Chrm.tours
        two_opt_on_route(tour, T, n_nodes)
    end
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function Find_Closeness(TT::Matrix{Float64}, h::Float64)
    n_nodes = size(TT)[1] - 2
    num = Int(ceil(h * n_nodes))
    ClosenessT = zeros(Int, n_nodes+1, num)
    @inbounds for i = 2:n_nodes+2
        a = copy(TT[i, 2:n_nodes+1])
        b = sortperm(a)
        ClosenessT[i-1, :] = b[2:num+1]
    end
    return ClosenessT
end