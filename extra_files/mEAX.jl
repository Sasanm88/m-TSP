using OffsetArrays

mutable struct abEdge
    edge::Tuple{Int, Int}
    first::Bool
end

mutable struct ABcycle
    edges::Vector{abEdge}
    length::Int
    nodes::Vector{Int}
end

function find_AB_cycles_m(p1::Vector{Tour}, p2::Vector{Tour}, n::Int)
    EA = Vector{Tuple{Int, Int}}()
    EAp = Vector{Tuple{Int, Int}}()
    EB = Vector{Tuple{Int, Int}}()

    for tour in p1
        temp1 = copy(tour.Sequence)
        
        push!(temp1, 0)
        pushfirst!(temp1, 0)
        EA = vcat(EA, [(temp1[i], temp1[i+1]) for i in 1:length(temp1)-1])
        EAp = vcat(EAp, [(temp1[i+1], temp1[i]) for i in 1:length(temp1)-1])
    end
    for tour in p2
        temp2 = copy(tour.Sequence)
        push!(temp2, 0)
        pushfirst!(temp2, 0)
        EB = vcat(EB, [(temp2[i], temp2[i+1]) for i in 1:length(temp2)-1])

    end
    ab = intersect(EA, EB)
    apb = intersect(EAp, EB)
    abp = Vector{Tuple{Int, Int}}()
    for pair in apb
        push!(abp, (pair[2], pair[1]))
    end
    EA_ = copy(EA)
    EA = setdiff(EA, ab, abp)
    EB = setdiff(EB, ab, apb)
#     return EA, EB
    paths = Vector{Vector{Int}}()
    AB_Cycles = ABcycle[]
    if isempty(EA)
        return AB_Cycles, EA
    end
    vertex = rand(collect(Iterators.flatten(EA)))
    first_parent = true
    path = [vertex]
    cycle = abEdge[]
    masks = zeros(Int, n+1)
    masks[vertex+1] = 1
    count = 1
    
    while !isempty(EA) || !isempty(EB)
        count += 1
        next_edge = (0,0)

        if first_parent
            if 1==2 #count%2 == 0
                linked_edges = filter(e -> (e[1] == vertex && e[2] != path[1]) || (e[2] == vertex && e[1] != path[1]), EA)
            else
                linked_edges = filter(e -> e[1] == vertex || e[2] == vertex, EA)
            end
            next_edge = rand(linked_edges)
            deleteat!(EA, findfirst(x->x==next_edge, EA))
        else
            if 1==2 #count%2 == 0
                linked_edges = filter(e -> (e[1] == vertex && e[2] != path[1]) || (e[2] == vertex && e[1] != path[1]), EB)
            else
                linked_edges = filter(e -> e[1] == vertex || e[2] == vertex, EB)
            end
            next_edge = rand(linked_edges)
            deleteat!(EB, findfirst(x->x==next_edge, EB))
        end

        if vertex == next_edge[2]
            vertex = next_edge[1]
        else
            vertex = next_edge[2]
        end
        push!(path, vertex)
        push!(cycle, abEdge(next_edge, first_parent))

#         println("edge: $next_edge from parent $(2-Int(first_parent)), vertex: $vertex, count: $count, masks[$(vertex)]=$(masks[vertex+1])")
        if masks[vertex+1] > 0 
            if (count - masks[vertex+1]) % 2 == 0
                temp = copy(cycle[masks[vertex+1]:count-1])
                temp2 = copy(path[masks[vertex+1]:count])
#                 for i in temp2
#                     print(i, " ")
#                 end
#                 println()
                push!(AB_Cycles, ABcycle(temp, count-masks[vertex+1], temp2[1:end-1]))
                push!(paths, temp2)
                
                if vertex == path[1] 
                    path = Int[]
                    cycle = abEdge[]
                    masks = zeros(Int, n+1)
#                     println("One complete tour, started and ended at $vertex")
                else
                    deleteat!(cycle, [i for i=masks[vertex+1]:count-1])
                    deleteat!(path, [i for i=masks[vertex+1]+1:count])
                    count = masks[vertex+1]
                    masks = zeros(Int, n+1)
                    for (i, node) in enumerate(path)
                        if masks[node+1] == 0 
                            masks[node+1] = i
                        end
                    end
                end

                if length(path) == 0
                    if first_parent
                        remaining_nodes = collect(Iterators.flatten(EB))
                    else
                        remaining_nodes = collect(Iterators.flatten(EA))
                    end
                    if isempty(remaining_nodes)
                        break
                    end
                    vertex = rand(remaining_nodes)
                    path = [vertex]
#                     println("Started over from node $vertex")
                    count = 1
                    masks[vertex+1] = count
                end
            end
        else
            masks[vertex+1] = count
        end
    #     if cycle[1] == vertex
    #         push!(cycles, cycle)
    #     end
        first_parent = !first_parent
    end
    return AB_Cycles , EA_
end


function separate_subtours_m(solutions::Vector{Tuple{Int, Int}}, n::Int, m::Int)
    neighbors = Dict{Int, Vector{Int}}()
    for edge in solutions
        if edge[1] in keys(neighbors)
            push!(neighbors[edge[1]], edge[2])
        else
            neighbors[edge[1]] = [edge[2]]
        end

        if edge[2] in keys(neighbors)
            push!(neighbors[edge[2]], edge[1])
        else
            neighbors[edge[2]] = [edge[1]]
        end
    end

    # Define a function to perform depth-first search
    function dfs1_m(vertex::Int, prev::Int, neighbors::Dict{Int64, Vector{Int64}}, visited::Vector{Bool}, tour::Vector{Int})
        visited[vertex] = true
        for neighbor in neighbors[vertex]
            if neighbor != prev
                if neighbor > 0
                    visited[neighbor] = true
                    push!(tour, neighbor)
                    return neighbor
                else
                    return neighbor
                end
            end
        end
    end

    function dfs_m(vertex::Int, neighbors::Dict{Int64, Vector{Int64}}, visited::Vector{Bool}, tour::Vector{Int})
        visited[vertex] = true
        for neighbor in neighbors[vertex]
            if !visited[neighbor]
                dfs_m(neighbor, neighbors, visited, tour)
            end
        end
        push!(tour, vertex)
    end

    # Find all subtours using depth-first search
    subtours = Vector{Vector{Int}}()
    visited = fill(false, n)
    for i=1:n
        if sum(neighbors[i]) == 0
            m = m-1
            push!(subtours, Int[0,i])
            visited[i] = true
            deleteat!(neighbors[0], findall(x->x==i, neighbors[0]))
        end
    end
    for i=1:m
        vertex = pop!(neighbors[0])
        prev = 0
        tour = Int[0, vertex]
        while true
            temp = dfs1_m(vertex, prev, neighbors, visited, tour)
            prev = vertex
            vertex = temp
            if temp == 0
                deleteat!(neighbors[0], findfirst(x->x==prev, neighbors[0]))
                break
            end
        end
        push!(subtours, tour)
    end
    remaining_nodes = findall(x->x==false, visited)
    for vertex in remaining_nodes
        if !visited[vertex]
            tour = Int[]
            dfs_m(vertex, neighbors, visited, tour)
            push!(subtours, reverse(tour))
        end
        visited[1] = false
    end
    return subtours
end


function update_cost_after_adding_the_subtour(tour::Tour, tour1::Vector{Int}, T::Matrix{Float64}, n1::Int, n2::Int, n3::Int, n4::Int)
    new_cost = sum(T[tour1[i]+1, tour1[i+1]+1] for i=1:length(tour1)-1) + T[tour1[1]+1, tour1[end]+1]
    tour.cost += new_cost
    tour.cost += -T[tour1[n1]+1, tour1[n3]+1] - T[tour.Sequence[n2]+1, tour.Sequence[n4]+1] + T[tour1[n1]+1, tour.Sequence[n2]+1] + T[tour1[n3]+1, tour.Sequence[n4]+1]
end

function modify_intermediate_solution(subtours::Vector{Vector{Int}}, T::Matrix{Float64})

    real_tours = Tour[]
    dlt_idx = Int[]
    for (i,tour) in enumerate(subtours)
        if tour[1] == 0
            push!(dlt_idx, i)
            push!(real_tours, Tour(tour, find_tour_length(tour[2:end], T)))
        end
    end

    deleteat!(subtours, dlt_idx)
    while length(subtours) > 0
        costs = [x.cost for x in real_tours]
        penalties = costs/maximum(costs)
#         penalties = ones(length(real_tours))
        sort!(subtours, by=x->length(x), rev=true)
        smallest_tour = pop!(subtours);
        best_n1 = 0
        best_n2 = 0 
        best_n3 = 0
        best_n4 = 0
        min_cost = Inf
        best_tour = 0
        if length(smallest_tour) == 2
            n1 = 1
            n3 = 2
            for i in 1:length(real_tours)
                tour = real_tours[i].Sequence
                if length(tour) == 2
                    n2 = 1
                    n4 = 2
                    new_cost = (T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                    if new_cost < min_cost
                        min_cost =  new_cost
                        best_tour = i
                        best_n1 = n1
                        best_n2 = n2 
                        best_n3 = n3
                        best_n4 = n4
                    end
                    n4 = 1
                    n2 = 2
                    new_cost = (T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                    if new_cost < min_cost
                        min_cost =  new_cost
                        best_tour = i
                        best_n1 = n1
                        best_n2 = n2 
                        best_n3 = n3
                        best_n4 = n4
                    end
                else
                    for j=1:length(tour)                    
                        if j == length(tour)
                            n2 = j
                            n4 = 1
                        else
                            n2 = j
                            n4 = j+1
                        end
                        cost = - T[tour[n2]+1, tour[n4]+1]
                        new_cost = (cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                        if new_cost < min_cost
                            min_cost =  new_cost
                            best_tour = i
                            best_n1 = n1
                            best_n2 = n2 
                            best_n3 = n3
                            best_n4 = n4
                        end
                        if j == length(tour)
                            n4 = j
                            n2 = 1
                        else
                            n2 = j+1
                            n4 = j
                        end
                        new_cost = (cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                        if new_cost < min_cost
                            min_cost =  new_cost
                            best_tour = i
                            best_n1 = n1
                            best_n2 = n2 
                            best_n3 = n3
                            best_n4 = n4
                        end
                    end
                end
            end
        else
            for n1=1:length(smallest_tour)
                if n1 == length(smallest_tour)
                    n3 = 1
                else
                    n3 = n1+1
                end
                cost1 = - T[smallest_tour[n1]+1, smallest_tour[n3]+1]
                for i in 1:length(real_tours)
                    tour = real_tours[i].Sequence
                    for j=1:length(tour)                    
                        if j == length(tour)
                            n2 = j
                            n4 = 1
                        else
                            n2 = j
                            n4 = j+1
                        end
                        cost = cost1 - T[tour[n2]+1, tour[n4]+1]
                        new_cost = (cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                        if new_cost < min_cost
                            min_cost =  new_cost
                            best_tour = i
                            best_n1 = n1
                            best_n2 = n2 
                            best_n3 = n3
                            best_n4 = n4
                        end
                        if j == length(tour)
                            n4 = j
                            n2 = 1
                        else
                            n2 = j+1
                            n4 = j
                        end
                        new_cost = (cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]) * penalties[i]
                        if new_cost < min_cost
                            min_cost =  new_cost
                            best_tour = i
                            best_n1 = n1
                            best_n2 = n2 
                            best_n3 = n3
                            best_n4 = n4
                        end
                    end
                end
            end
        end

        if best_n3 == 1
            if best_n4 == 1 && best_n2 == length(real_tours[best_tour].Sequence)
                new_tour = vcat(smallest_tour, reverse(real_tours[best_tour].Sequence))
            elseif best_n2 == 1 && best_n4 == length(real_tours[best_tour].Sequence)
                new_tour = vcat(smallest_tour, real_tours[best_tour].Sequence)
            else
                if best_n2 < best_n4
                    new_tour = vcat(reverse(smallest_tour), real_tours[best_tour].Sequence[best_n4:end], real_tours[best_tour].Sequence[1:best_n2])
                else
                    new_tour = vcat(smallest_tour, real_tours[best_tour].Sequence[best_n2:end], real_tours[best_tour].Sequence[1:best_n4])
                end
            end
        else
            if best_n4 == 1 && best_n2 == length(real_tours[best_tour].Sequence)
                new_tour = vcat(smallest_tour[1:best_n1], reverse(real_tours[best_tour].Sequence), smallest_tour[best_n3:end])
            elseif best_n2 == 1 && best_n4 == length(real_tours[best_tour].Sequence)
                new_tour = vcat(smallest_tour[1:best_n1], real_tours[best_tour].Sequence, smallest_tour[best_n3:end])
            else
                if best_n2 < best_n4
                    new_tour = vcat(smallest_tour[1:best_n1], reverse(real_tours[best_tour].Sequence[1:best_n2]), reverse(real_tours[best_tour].Sequence[best_n4:end]), smallest_tour[best_n3:end])
                else
                    new_tour = vcat(smallest_tour[1:best_n1], real_tours[best_tour].Sequence[best_n2:end], real_tours[best_tour].Sequence[1:best_n4], smallest_tour[best_n3:end])
                end
            end
        end
        update_cost_after_adding_the_subtour(real_tours[best_tour], smallest_tour, T, best_n1, best_n2, best_n3, best_n4)
        real_tours[best_tour].Sequence = new_tour
    end


    for tour in real_tours
        i = findfirst(x->x==0, tour.Sequence)
        if i == 1 
            popfirst!(tour.Sequence)
        elseif i == length(tour.Sequence)
            pop!(tour.Sequence)
        else
            tour.Sequence = vcat(tour.Sequence[i+1:end], tour.Sequence[1:i-1])
        end
    end
    return real_tours
end

function m_eax_crossover(p1::Vector{Tour}, p2::Vector{Tour}, T::Matrix{Float64}, n::Int)
    m = length(p1)
    abcycles, a = find_AB_cycles_m(p1, p2, n)
    if isempty(abcycles)
        return Int[]
    end
    select_effective_cycles(abcycles)
    
    cycle1 = popfirst!(abcycles)
    Eset = copy(cycle1.edges)
#     r = rand(1:length(abcycles))
#     cycle1 = abcycles[r]
#     deleteat!(abcycles, r)
#     Eset = copy(cycle1.edges)
    existing_nodes = copy(cycle1.nodes)
    for cycle in abcycles
        if any(x -> x in existing_nodes, cycle.nodes)
            Eset = vcat(Eset, cycle.edges)
            existing_nodes = union(existing_nodes, cycle.nodes)
        end
    end

    solutions = generate_intermediate_solution(a, Eset)
    subtours = separate_subtours_m(solutions, n, m)

    tours = modify_intermediate_solution(subtours, T)
    chrm = Chromosome(Int[], 0.0, 0.0, tours)
    for tour in tours
        chrm.genes = vcat(chrm.genes, tour.Sequence)
        if tour.cost > chrm.fitness
            chrm.fitness = tour.cost
        end
    end
    return chrm.genes
end
