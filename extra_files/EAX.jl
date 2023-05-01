mutable struct abEdge
    edge::Tuple{Int, Int}
    first::Bool
end

mutable struct ABcycle
    edges::Vector{abEdge}
    length::Int
end

function find_AB_cycles(p1::Vector{Int}, p2::Vector{Int})
    temp1 = copy(p1)
    temp2 = copy(p2)
    push!(temp1, 0)
    pushfirst!(temp1, 0)
    push!(temp2, 0)
    pushfirst!(temp2, 0)

    EA_ = [(temp1[i], temp1[i+1]) for i in 1:length(temp1)-1]
    EB = [(temp2[i], temp2[i+1]) for i in 1:length(temp2)-1];
    EA = copy(EA_)
    paths = Vector{Vector{Int}}()
    AB_Cycles = ABcycle[]
    vertex = rand(0:n)
    first_parent = true
    path = [vertex]
    cycle = abEdge[]
    masks = zeros(Int, n+1)
    masks[vertex+1] = 1
    count = 1
    println(vertex)
    while !isempty(EA) || !isempty(EB)
        count += 1
        next_edge = (0,0)
        while true
            if first_parent
                linked_edges = filter(e -> e[1] == vertex || e[2] == vertex, EA) #if empty select from other parent
                if isempty(linked_edges)
                    first_parent = !first_parent
                else
                    next_edge = rand(linked_edges)
                    deleteat!(EA, findfirst(x->x==next_edge, EA))
                    break
                end
            else
                linked_edges = filter(e -> e[1] == vertex || e[2] == vertex, EB)
                if isempty(linked_edges)
                    first_parent = !first_parent
                else
                    next_edge = rand(linked_edges)
                    deleteat!(EB, findfirst(x->x==next_edge, EB))
                    break
                end
            end
        end

        if vertex == next_edge[2]
            vertex = next_edge[1]
        else
            vertex = next_edge[2]
        end
        push!(path, vertex)
        push!(cycle, abEdge(next_edge, first_parent))

    #     println("edge: $next_edge, vertex: $vertex, count: $count, masks[$(vertex)]=$(masks[vertex+1])")
        if masks[vertex+1] > 0 
            temp = copy(cycle[masks[vertex+1]:count-1])
            temp2 = copy(path[masks[vertex+1]:count])
            push!(AB_Cycles, ABcycle(temp, count-masks[vertex+1]))
            push!(paths, temp2)
            if vertex == path[1] 
                path = Int[]
                cycle = abEdge[]
                masks = zeros(Int, n+1)
    #             println("One complete tour, started and ended at $vertex")
            else
                deleteat!(cycle, [i for i=masks[vertex+1]:count-1])
                deleteat!(path, [i for i=masks[vertex+1]+1:count])
                count = masks[vertex+1]
                for i in temp2[2:length(temp2)]
                    masks[i+1] = 0
                end
            end

            if length(path) == 0

                remaining_nodes = collect(intersect(Iterators.flatten(EA), Iterators.flatten(EB)))
                if isempty(remaining_nodes)
                    break
                end
                vertex = rand(remaining_nodes)
                path = [vertex]
    #             println("Started over from node $vertex")
                count = 1
                masks[vertex+1] = count
            end
        else
            masks[vertex+1] = count
        end
    #     if cycle[1] == vertex
    #         push!(cycles, cycle)
    #     end
        first_parent = !first_parent
    end
    return AB_Cycles, EA_
end

function select_effective_cycles(AB_Cycles::Vector{ABcycle})
    delete_idx = Int[]
    for (i,cycle) in enumerate(AB_Cycles)
        if cycle.length == 2
            push!(delete_idx, i)
        end
    end
    deleteat!(AB_Cycles, delete_idx)
    sort!(AB_Cycles, by=x->x.length, rev=true)
end
    

function form_Eset_rand(AB_Cycles::Vector{ABcycle})
    Eset = abEdge[]
    for cycle in AB_Cycles
        if rand() < 0.5
            Eset = vcat(Eset, cycle.edges)
        end
    end
    if isempty(Eset)
        Eset = rand(AB_Cycles).edges
    end
    return Eset
end

function generate_intermediate_solution(tourA::Vector{Tuple{Int, Int}}, Eset::Vector{abEdge})
    Aset = Vector{Tuple{Int, Int}}()
    Bset = Vector{Tuple{Int, Int}}()
    for ab in Eset
        if ab.first
            push!(Aset, ab.edge)
        else
            push!(Bset, ab.edge)
        end
    end
    temp = setdiff(tourA, Aset)
    return vcat(temp, Bset)
end

function separate_subtours(solutions::Vector{Tuple{Int, Int}})
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
    function dfs(vertex::Int, neighbors::Dict{Int64, Vector{Int64}}, visited::Vector{Bool}, tour::Vector{Int})
        visited[vertex+1] = true
        for neighbor in neighbors[vertex]
            if !visited[neighbor+1]
                dfs(neighbor, neighbors, visited, tour)
            end
        end
        push!(tour, vertex)
    end

    # Find all subtours using depth-first search
    subtours = Vector{Vector{Int}}()
    visited = fill(false, n+1)
    for vertex in keys(neighbors)
        if !visited[vertex+1]
            tour = Int[]
            dfs(vertex, neighbors, visited, tour)
            push!(subtours, reverse(tour))
        end
    end
    return subtours
end

function modify_intermediate_solution(subtours::Vector{Vector{Int}}, T::Matrix{Float64})
    while length(subtours) > 1
        sort!(subtours, by=x->length(x), rev=true)
        smallest_tour = pop!(subtours);
        best_n1 = 0
        best_n2 = 0 
        best_n3 = 0
        best_n4 = 0
        min_cost = Inf
        best_tour = 0
        for n1=1:length(smallest_tour)-1
            n3 = n1+1
            cost = - T[smallest_tour[n1]+1, smallest_tour[n3]+1]
            for (i,tour) in enumerate(subtours)
                for j=1:length(tour)-1
                    cost -= T[tour[j]+1, tour[j+1]+1]
                    n2 = j
                    n4 = j+1
                    new_cost = cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
                    if new_cost < min_cost
                        min_cost =  new_cost
                        best_tour = i
                        best_n1 = n1
                        best_n2 = n2 
                        best_n3 = n3
                        best_n4 = n4
                    end
                    n2 = j+1
                    n4 = j
                    new_cost = cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
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

        if best_n2 < best_n4
            new_tour = vcat(smallest_tour[1:best_n1], reverse(subtours[best_tour][1:best_n2]), reverse(subtours[best_tour][best_n4:end]), smallest_tour[best_n3:end])
        else
            new_tour = vcat(smallest_tour[1:best_n1], subtours[best_tour][best_n2:end], subtours[best_tour][1:best_n4], smallest_tour[best_n3:end])
        end
        subtours[best_tour] = new_tour
    end
    new_tour = subtours[1]
    i = findfirst(x->x==0, new_tour)
    if i == 1 
        return new_tour[2:end]
    elseif i == length(new_tour)
        return new_tour[1:end-1]
    else
        return vcat(new_tour[i+1:end], new_tour[1:i-1])
    end
end

function EAX_rand(p1::Vector{Int}, p2::Vector{Int})
    AB_Cycles, EA = find_AB_cycles(p1, p2)
    select_effective_cycles(AB_Cycles)
    Eset = form_Eset_rand(AB_Cycles)
    solutions = generate_intermediate_solution(EA, Eset)
    subtours = separate_subtours(solutions)
    new_tour = modify_intermediate_solution(subtours, T)
    return new_tour
end