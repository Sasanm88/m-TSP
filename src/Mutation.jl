function two_opt_mutation(Chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int)   #2-opt 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1, k2 = sort(sample(1:length(tour1), 2, replace = false))

    new_cost = Calculate_new_cost_2_opt(tour1, cost1, k1, k2, T, n_nodes)

    tour1[k1:k2] = reverse(tour1[k1:k2])
    Chrm.tours[r1].cost = new_cost
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function cross_mutation(Chrm::Chromosome, Customers::Matrix{Float64}, depot::Vector{Float64}, T::Matrix{Float64}, n_nodes::Int)
    m = length(Chrm.tours)
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour_neighbors = find_tour_neighbors(Chrm.tours, Customers, depot, m)
    r2 = tour_neighbors[r1][rand(1:length(tour_neighbors[r1]))]
    
    t1 = Chrm.tours[r1].Sequence
    t2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    n1 = length(t1)
    n2 = length(t2)
    
    if n1 < 2 || n2 < 2
        return Chrm
    end
    
    k11, k12 = sort(sample(1:n1, 2, replace = false))
    k21, k22 = sort(sample(1:n2, 2, replace = false))
    
    new_cost1, new_cost2, straight1, straight2 = Calculate_new_cost_cross(t1, cost1, t2, cost2, k11, k12, k21, k22, T, n_nodes)

    if straight2
        alpha1 = copy(t1[k11:k12])
    else
        alpha1 = reverse(copy(t1[k11:k12]))
    end
    if straight1
        alpha2 = copy(t2[k21:k22])
    else
        alpha2 = reverse(copy(t2[k21:k22]))
    end
    deleteat!(t1, [i for i=k11:k12])
    for i=1:k22-k21+1
        insert!(t1, i+k11-1, alpha2[i])
    end

    deleteat!(t2, [i for i=k21:k22])
    for i=1:k12-k11+1
        insert!(t2, i+k21-1, alpha1[i])
    end

    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function mix_neighbors_mutation(Chrm::Chromosome, Customers::Matrix{Float64}, depot::Vector{Float64}, T::Matrix{Float64}, n_nodes::Int)
    m = length(Chrm.tours)
    
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 4
        return Chrm
    end
    tour_neighbors = find_tour_neighbors(Chrm.tours, Customers, depot, m)
    r2 = tour_neighbors[r1][1]
    tour2 = Chrm.tours[r2].Sequence
    if length(tour2) <= 4
        return Chrm
    end
    
    if length(tour1) <= length(tour2)
        idx1 , idx2 = sort(sample(2:length(tour1)-1, 2, replace = false))
        t1 = vcat(tour2[1:idx1-1], tour1[idx1:idx2], tour2[idx2+1:length(tour2)])
        t2 = vcat(tour1[1:idx1-1], tour2[idx1:idx2], tour1[idx2+1:length(tour1)])
    else
        idx1 , idx2 = sort(sample(2:length(tour2)-1, 2, replace = false))
        t1 = vcat(tour1[1:idx1-1], tour2[idx1:idx2], tour1[idx2+1:length(tour1)])
        t2 = vcat(tour2[1:idx1-1], tour1[idx1:idx2], tour2[idx2+1:length(tour2)])           
    end    
    Chrm.tours[r1] = Tour(t1, find_tour_length(t1, T))
    Chrm.tours[r2] = Tour(t2, find_tour_length(t2, T))
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end



function scatter_mutation(Chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int)
    m = length(Chrm.tours)
    
    moving_nodes = Int[]
  
    for tour in Chrm.tours
        c = sort(sample(1:length(tour.Sequence), rand(1:min(3,length(tour.Sequence))), replace = false))
        cities = copy(tour.Sequence[c])
#         Remove_cities_from_one_tour(tour, c, T, n_nodes)
        deleteat!(tour.Sequence, c)
        tour.cost = find_tour_length(tour.Sequence, T)
        moving_nodes = vcat(moving_nodes, cities)
    end
    for city in moving_nodes 
        put_city_in_tour(Chrm.tours, city, T, n_nodes)
    end
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function mutate(Chrm::Chromosome, Customers::Matrix{Float64}, depot::Vector{Float64}, T::Matrix{Float64}, n_nodes::Int)
    
    r = rand()
    if r < 0.3
        return two_opt_mutation(Chrm, T, n_nodes)
    elseif r < 0.7
        return cross_mutation(Chrm, Customers, depot, T, n_nodes)
    else
        return scatter_mutation(Chrm, T, n_nodes)
    end
        
end



# function Displacement_mutation(child::Vector{Int}, n_nodes::Int)  
#     idx1 = rand(2:n_nodes-1)
#     idx2 = rand(2:n_nodes-1)
#     while idx1 == idx2
#         idx1 = rand(2:n_nodes-1)
#         idx2 = rand(2:n_nodes-1)
#     end

#     if idx1 > idx2
#         temp = idx1
#         idx1 = idx2
#         idx2 = temp
#     end

#     subtour = child[idx1:idx2]
#     deleteat!(child, idx1:idx2)
#     idx3 = rand(1:length(child))

#     for i = 1:idx2-idx1+1
#         insert!(child, idx3 + i, subtour[i])
#     end
# end

# function Exchange_mutation(child::Vector{Int}, n_nodes::Int)
#     idx1, idx2 = sample(1:n_nodes, Weights(ones(n_nodes)), 2, replace=false)
#     child[idx1], child[idx2] = child[idx2], child[idx1]
# end


# function Insertion_mutation(child::Vector{Int}, n_nodes::Int)
#     idx1 = rand(1:n_nodes)
#     city = child[idx1]
#     deleteat!(child, idx1)
#     idx2 = rand(1:length(child))
#     insert!(child, idx2, city)
# end

# function Simple_Inversion_mutation(child::Vector{Int}, n_nodes::Int)
#     idx1 = rand(1:n_nodes)
#     idx2 = rand(1:n_nodes)
#     while idx1 == idx2
#         idx1 = rand(1:n_nodes)
#         idx2 = rand(1:n_nodes)
#     end

#     if idx1 > idx2
#         temp = idx1
#         idx1 = idx2
#         idx2 = temp
#     end
#     child[idx1:idx2] = reverse(child[idx1:idx2])
# end


# function Inversion_mutation(child::Vector{Int}, n_nodes::Int)
#     idx1 = rand(2:n_nodes-1)
#     idx2 = rand(2:n_nodes-1)
#     while idx1 == idx2
#         idx1 = rand(2:n_nodes-1)
#         idx2 = rand(2:n_nodes-1)
#     end

#     if idx1 > idx2
#         temp = idx1
#         idx1 = idx2
#         idx2 = temp
#     end

#     subtour = child[idx1:idx2]
#     deleteat!(child, idx1:idx2)
#     idx3 = rand(1:length(child))
#     reverse!(subtour)
#     for i = 1:idx2-idx1+1
#         insert!(child, idx3 + i, subtour[i])
#     end
# end


# function Scramble_mutation(child::Vector{Int}, n_nodes::Int)
#     idx1 = rand(1:n_nodes)
#     idx2 = rand(1:n_nodes)
#     while idx1 == idx2
#         idx1 = rand(1:n_nodes)
#         idx2 = rand(1:n_nodes)
#     end

#     if idx1 > idx2
#         temp = idx1
#         idx1 = idx2
#         idx2 = temp
#     end
#     child[idx1:idx2] = shuffle!(child[idx1:idx2])
# end

# function tour_mutation(child::Vector{Int}, n_nodes::Int)
#     idx = sample(1:n_nodes, Int(round(0.2 * n_nodes)), replace=false)
#     child[idx] = child[shuffle(idx)]
# end

function find_tour_length(tt::Vector{Int}, T::Matrix{Float64})
    t = copy(tt)
    pushfirst!(t, 0)
    push!(t, 0)
    z = 0.0
    for i = 1:length(t)-1
        z += T[t[i]+1, t[i+1]+1]
    end
    return z
end

# function new_mutation(child::Chromosome, T::Matrix{Float64}, pm::Float64)
#     removed_cities = Int[]
#     for tour in child.tours
#         indices = Int[]
#         for i = 1:length(tour.Sequence)
#             if rand() < pm
#                 push!(indices, i)
#                 push!(removed_cities, tour.Sequence[i])
#             end
#         end
#         deleteat!(tour.Sequence, indices)
#         tour.cost = find_tour_length(tour.Sequence, T)
#     end
#     for city in removed_cities
#         put_city_in_tour(child.tours, city, T, n_nodes)
#     end
#     child.fitness = maximum([tour.cost for tour in child.tours])
#     return child
# end

# function Mutate(child::Vector{Int}, Mutation_Chance::Float64)
#     n_nodes = length(child)
#     Mutation_functions = [Displacement_mutation, Simple_Inversion_mutation, Inversion_mutation]
#     if rand() < Mutation_Chance
#         Mutation_functions[rand(1:length(Mutation_functions))](child, n_nodes)
#     end
# end

# function chunk_mutation(chrm::Chromosome, T::Matrix{Float64}, chunk_length::Int, n_nodes::Int)
#     m = length(chrm.tours)
#     chunks = Vector{Vector{Int}}()
#     for tour in chrm.tours
#         nt = length(tour.Sequence)
#         for i = 1:Int(ceil(nt/chunk_length))
#             push!(chunks, tour.Sequence[(i-1)*chunk_length+1:min(nt,i*chunk_length)])
#         end
#     end
#     new_tours = Tour[]
#     for i=1:m
#         push!(new_tours, Tour(Int[],0.0))
#     end 

#     last_added_tour = zeros(m) 

#     while length(chunks) > 0
#         best_chunk = 0
#         least_distance = Inf
#         best_tour = 0
#         reverse_ = false
#         tours_indices = findall(x->x==0, last_added_tour)
#         for i in tours_indices
#             nt = length(new_tours[i].Sequence)
#             last_node = 0
#             if nt > 0
#                 last_node = new_tours[i].Sequence[nt]
#             end
#             for (j,chunk) in enumerate(chunks)
#                 if length(chunks) > m
#                     if T[last_node+1, chunk[1]+1] < least_distance
#                         best_tour = i
#                         best_chunk = j
#                         least_distance = T[last_node+1, chunk[1]+1]
#                         reverse_ = false
#                     elseif T[last_node+1, chunk[length(chunk)]+1] < least_distance
#                         best_tour = i
#                         best_chunk = j
#                         least_distance = T[last_node+1, chunk[length(chunk)]+1]
#                         reverse_ = true
#                     end
#                 else
#                     if T[last_node+1, chunk[1]+1] + T[chunk[length(chunk)]+1, n_nodes+2] < least_distance
#                         best_tour = i
#                         best_chunk = j
#                         least_distance = T[last_node+1, chunk[1]+1] + T[chunk[length(chunk)]+1, n_nodes+2]
#                         reverse_ = false
#                     elseif T[last_node+1, chunk[length(chunk)]+1] + T[chunk[1]+1, n_nodes+2] < least_distance
#                         best_tour = i
#                         best_chunk = j
#                         least_distance = T[last_node+1, chunk[length(chunk)]+1] + T[chunk[1]+1, n_nodes+2]
#                         reverse_ = true
#                     end
#                 end
#             end
#         end
#         if reverse_
#             new_tours[best_tour].Sequence = vcat(new_tours[best_tour].Sequence, reverse(chunks[best_chunk]))
#         else
#             new_tours[best_tour].Sequence = vcat(new_tours[best_tour].Sequence, chunks[best_chunk])
#         end
#         last_added_tour[best_tour] += 1
#         deleteat!(chunks, best_chunk)
#         if sum(last_added_tour) == m
#             last_added_tour = zeros(m)
#         end
#     end
#     new_chrm = Chromosome(Int[], 0.0, 0.0, new_tours)
#     for tour in new_tours
# #         t1 = copy(tour.Sequence)
# #         pushfirst!(t1, 0)
# #         t2, z2 = find_tsp_tour1(T[t1.+1, t1.+1])
#         tour.cost = find_tour_length(tour.Sequence, T)
# #         tour.cost = z2
# #         tour.Sequence = tour.Sequence[t2]
#         new_chrm.genes = vcat(new_chrm.genes, tour.Sequence)
#         if tour.cost > new_chrm.fitness
#             new_chrm.fitness = tour.cost
#         end
#     end
#     return new_chrm
# end

# function chunk_mutation_rand(chrm::Chromosome, T::Matrix{Float64}, chunk_length::Int)
#     m = length(chrm.tours)
#     chunks = Vector{Vector{Int}}()
#     for tour in chrm.tours
#         nt = length(tour.Sequence)
#         for i = 1:Int(ceil(nt/chunk_length))
#             push!(chunks, tour.Sequence[(i-1)*chunk_length+1:min(nt,i*chunk_length)])
#         end
#     end
#     new_tours = Tour[]
#     for i=1:m
#         push!(new_tours, Tour(Int[],0.0))
#     end 

#     last_added_tour = zeros(m) 

#     while length(chunks) > 0
#         best_chunk = rand(1:length(chunks))
#         reverse_ = true
#         if rand() < 0.5
#             reverse_ = false
#         end
#         tours_indices = findall(x->x==0, last_added_tour)
#         best_tour = tours_indices[rand(1:length(tours_indices))]
#         if reverse_
#             new_tours[best_tour].Sequence = vcat(new_tours[best_tour].Sequence, reverse(chunks[best_chunk]))
#         else
#             new_tours[best_tour].Sequence = vcat(new_tours[best_tour].Sequence, chunks[best_chunk])
#         end
#         last_added_tour[best_tour] += 1
#         deleteat!(chunks, best_chunk)
#         if sum(last_added_tour) == m
#             last_added_tour = zeros(m)
#         end
#     end
#     new_chrm = Chromosome(Int[], 0.0, 0.0, new_tours)
#     for tour in new_tours
#         t1 = copy(tour.Sequence)
#         pushfirst!(t1, 0)
#         t2, z2 = find_tsp_tour1(T[t1.+1, t1.+1])
# #         tour.cost = find_tour_length(tour.Sequence, T)
#         tour.cost = z2
#         tour.Sequence = tour.Sequence[t2]
#         new_chrm.genes = vcat(new_chrm.genes, tour.Sequence)
#         if tour.cost > new_chrm.fitness
#             new_chrm.fitness = tour.cost
#         end
#     end
#     return new_chrm
# end

# # function put_one_city_in_one_tour(c::Tour, city::Int, T::Matrix{Float64}, n_nodes::Int)
# #     least_increase = Inf
# #     best_position = 0
# #     tour = c.Sequence
# #     nt = length(tour)
# #     if nt==0
# #         increase = T[1, city+1] + T[city+1, n_nodes+2]
# #         best_position = 1
# #     else
# #         increase = T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1]
# #         if increase < least_increase
# #             least_increase = increase
# #             best_position = 1
# #         end
# #         for j = 2:nt
# #             increase = T[tour[j-1]+1, city+1] + T[city+1, tour[j]+1] - T[tour[j-1]+1, tour[j]+1]
# #             if increase < least_increase
# #                 least_increase = increase
# #                 best_position = j
# #             end
# #         end
# #         increase = T[tour[nt]+1, city+1] + T[city+1, n_nodes+2] - T[tour[nt]+1, n_nodes+2]
# #         if increase < least_increase
# #             least_increase = increase
# #             best_position = nt+1
# #         end
# #     end
# #     insert!(tour, best_position, city)
# #     c.cost += least_increase
# # end 

# function put_city_in_one_tour(c::Tour, city::Int, T::Matrix{Float64}, n_nodes::Int)
#     least_increase = Inf
#     best_position = 0
#     tour = c.Sequence
#     nt = length(tour)
#     if nt==0
#         least_increase = T[1, city+1] + T[city+1, n_nodes+2]
#         best_position = 1
#     else
#         increase = T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1]
#         if increase < least_increase
#             least_increase = increase
#             best_position = 1
#         end
#         for j = 2:nt
#             increase = T[tour[j-1]+1, city+1] + T[city+1, tour[j]+1] - T[tour[j-1]+1, tour[j]+1]
#             if increase < least_increase
#                 least_increase = increase
#                 best_position = j
#             end
#         end
#         increase = T[tour[nt]+1, city+1] + T[city+1, n_nodes+2] - T[tour[nt]+1, n_nodes+2]
#         if increase < least_increase
#             least_increase = increase
#             best_position = nt+1
#         end
#     end
#     insert!(c.Sequence, best_position, city)
#     c.cost += least_increase
# end   

# function scatter_mutation(chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int, m::Int, demands::Vector{Int}, W::Int, Customers::Matrix{Float64}, depot::Vector{Float64}, h1::Float64)
#     t1 = chrm.tours[1].Sequence;
#     sort!(chrm.tours, by=x->x.cost, rev=true)
#     means = [mean(Customers[t1, :], dims=1)[1,:] for t1 in [chrm.tours[i].Sequence for i=1:length(chrm.tours)]]

#     distances = zeros(m, m)
#     for i = 1:m-1
#         for j = i+1:m
#             distances[i,j] = euclidean(means[i], means[j])
#             distances[j,i] = distances[i,j]
#         end
#     end

#     all_tours = [i for i=1:m]

#     current_tour = 1
#     while length(all_tours) > 0
#         next_tour = 0
#         least_dist = 100000

#         if length(all_tours) == 1
#             next_tour = 1
#         else
#             for i in all_tours
#                 if distances[current_tour, i] < least_dist
#                     if distances[current_tour, i] > 0
#                         least_dist = distances[current_tour, i]
#                         next_tour = i
#                     end
#                 end
#             end
#         end

#         t1 = copy(chrm.tours[current_tour].Sequence)
#         tau = rand(1:Int(round(h1*length(t1))))
        
#         Candidates = t1[sortperm([euclidean(Customers[i, :], means[next_tour]) for i in t1])][1:tau]

#         cities = findall(x->x in Candidates, t1)
#         Remove_cities_from_one_tour(chrm.tours[current_tour], cities, T, n_nodes)
#         deleteat!(chrm.tours[current_tour].Sequence, sort(cities))
#         for city in Candidates
#             put_city_in_one_tour(chrm.tours[next_tour], city, T, n_nodes)
#         end

#         deleteat!(all_tours, findfirst(x->x==current_tour, all_tours))
#         current_tour = next_tour
#     end
#     chrm.genes = Int[]
#     chrm.fitness = 0
#     for tour in chrm.tours
#         chrm.genes = vcat(chrm.genes, tour.Sequence)
#         if chrm.fitness < tour.cost
#             chrm.fitness = tour.cost
#         end
#     end
#     obj, trips = SPLIT(T, demands, m, W, chrm.genes)
#     chrm1 = Chromosome(chrm.genes, obj, 0.0, trips)
#     return chrm1
# end

# function prob_mutation(chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int, p0::Float64)
#     tour1 = chrm.tours[1]
#     tour2 = chrm.tours[2]
#     t1 = copy(tour1.Sequence)
#     t2 = copy(tour2.Sequence)
#     D = zeros(length(t1),length(t2))
#     for i=1:length(t1)
#         for j=1:length(t2)
#             D[i,j] = T[t1[i]+1, t2[j]+1]
#         end
#     end
#     p1 = zeros(length(t1))
#     p2 = zeros(length(t2))
    
#     for (i, node) in enumerate(sortperm(minimum(D, dims=2)[:,1]))
#         p1[node] = p0/i
#     end
#     for (i, node) in enumerate(sortperm(minimum(D, dims=1)[1,:]))
#         p2[node] = p0/i
#     end
#     delete1 = Int[]
#     delete2 = Int[]
#     for i=1:length(p1)
#         if rand() < p1[i]
#             push!(delete1, i)
#         end
#     end
#     for i=1:length(p2)
#         if rand() < p2[i]
#             push!(delete2, i)
#         end
#     end
#     if length(delete1) > 0
#         Remove_cities_from_one_tour(tour1, delete1, T, n_nodes)
#         deleteat!(tour1.Sequence, delete1)
#     end
#     if length(delete2) > 0
#         Remove_cities_from_one_tour(tour2, delete2, T, n_nodes)
#         deleteat!(tour2.Sequence, delete2)
#     end
#     for i in delete1
#         put_city_in_one_tour(tour2, t1[i], T, n_nodes)
#     end
#     for i in delete2
#         put_city_in_one_tour(tour1, t2[i], T, n_nodes)
#     end
#     chrm.genes = Int[]
#     chrm.fitness = maximum([chrm.tours[i].cost for i=1:length(chrm.tours)])
#     for tour in chrm.tours
#         chrm.genes = vcat(chrm.genes, tour.Sequence)
#     end
#     return chrm
# end


# function put_one_seq_in_one_tour(c::Tour, seq::Vector{Int}, T::Matrix{Float64}, n_nodes::Int)
#     least_increase = Inf
#     best_position = 0
#     tour = c.Sequence
#     nt = length(tour)
#     if nt==0
#         increase = T[1, city+1] + T[city+1, n_nodes+2]
#         best_position = 1
#     else
#         increase = T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1]
#         if increase < least_increase
#             least_increase = increase
#             best_position = 1
#         end
#         for j = 2:nt
#             increase = T[tour[j-1]+1, city+1] + T[city+1, tour[j]+1] - T[tour[j-1]+1, tour[j]+1]
#             if increase < least_increase
#                 least_increase = increase
#                 best_position = j
#             end
#         end
#         increase = T[tour[nt]+1, city+1] + T[city+1, n_nodes+2] - T[tour[nt]+1, n_nodes+2]
#         if increase < least_increase
#             least_increase = increase
#             best_position = nt+1
#         end
#     end
#     insert!(tour, best_position, city)
#     c.cost += least_increase
# end 

# function prob_chunk_mutation(chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int, p0::Float64)
#     tour1 = chrm.tours[1]
#     tour2 = chrm.tours[2]
#     t1 = copy(tour1.Sequence)
#     t2 = copy(tour2.Sequence)
#     D = zeros(length(t1),length(t2))
#     for i=1:length(t1)
#         for j=1:length(t2)
#             D[i,j] = T[t1[i]+1, t2[j]+1]
#         end
#     end
#     p1 = zeros(length(t1))
#     p2 = zeros(length(t2))
    
#     for (i, node) in enumerate(sortperm(minimum(D, dims=2)[:,1]))
#         p1[node] = p0/i
#     end
#     for (i, node) in enumerate(sortperm(minimum(D, dims=1)[1,:]))
#         p2[node] = p0/i
#     end
#     delete1 = Int[]
#     delete2 = Int[]
#     for i=1:length(p1)
#         if rand() < p1[i]
#             push!(delete1, i)
#         end
#     end
#     for i=1:length(p2)
#         if rand() < p2[i]
#             push!(delete2, i)
#         end
#     end
#     if length(delete1) > 0
#         Remove_cities_from_one_tour(tour1, delete1, T, n_nodes)
#         deleteat!(tour1.Sequence, delete1)
#     end
#     if length(delete2) > 0
#         Remove_cities_from_one_tour(tour2, delete2, T, n_nodes)
#         deleteat!(tour2.Sequence, delete2)
#     end
#     for i in delete1
#         put_city_in_one_tour(tour2, t1[i], T::Matrix{Float64}, n_nodes::Int)
#     end
#     for i in delete2
#         put_city_in_one_tour(tour1, t2[i], T::Matrix{Float64}, n_nodes::Int)
#     end
#     chrm.genes = Int[]
#     chrm.fitness = maximum([chrm.tours[i].cost for i=1:length(chrm.tours)])
#     for tour in chrm.tours
#         chrm.genes = vcat(chrm.genes, tour.Sequence)
#     end
#     return chrm
# end

# function mutation_cross(Chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int)   #Cross Exchange
#     r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
#     routes = [i for i=1:length(Chrm.tours)]
#     r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
#     t1 = Chrm.tours[r1].Sequence
#     t2 = Chrm.tours[r2].Sequence
#     cost1 = Chrm.tours[r1].cost
#     cost2 = Chrm.tours[r2].cost
#     if length(t1) < 6 || length(t2) < 6
#         return Chrm
#     end
#     tau = min(8, Int(round(min(length(t1), length(t2))/2)))
#     k11 = rand(1:length(t1)-tau)
#     l1 = rand(1:tau)
#     k12 = k11 + l1

#     k21 = rand(1:length(t2)-tau)
#     l2 = rand(1:tau)
#     k22 = k21 + l2
    
#     new_cost1, new_cost2, straight1, straight2 = Calculate_new_cost_cross(t1, cost1, t2, cost2, k11, k12, k21, k22, T, n_nodes)

#     if straight2
#         alpha1 = copy(t1[k11:k12])
#     else
#         alpha1 = reverse(copy(t1[k11:k12]))
#     end
#     if straight1
#         alpha2 = copy(t2[k21:k22])
#     else
#         alpha2 = reverse(copy(t2[k21:k22]))
#     end
#     deleteat!(t1, [i for i=k11:k12])
#     for i=1:k22-k21+1
#         insert!(t1, i+k11-1, alpha2[i])
#     end

#     deleteat!(t2, [i for i=k21:k22])
#     for i=1:k12-k11+1
#         insert!(t2, i+k21-1, alpha1[i])
#     end

#     Chrm.tours[r1].cost = new_cost1
#     Chrm.tours[r2].cost = new_cost2
#     Chrm.genes = Int[]
#     Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
#     for tour in Chrm.tours
#         Chrm.genes = vcat(Chrm.genes, tour.Sequence)
#     end
#     return Chrm
# end