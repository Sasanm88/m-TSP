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

function find_AB_cycles(p1::Vector{Int}, p2::Vector{Int}, n::Int)
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

function form_eset_1ab(AB_Cycles::Vector{ABcycle})
    Eset = AB_Cycles[1].edges
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

function separate_subtours(solutions::Vector{Tuple{Int, Int}}, n::Int)
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
        if length(smallest_tour) == 2
            n1 = 1
            n3 = 2
            for (i,tour) in enumerate(subtours)
                if length(tour) == 2
                    n2 = 1
                    n4 = 2
                    new_cost = T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
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
                    new_cost = T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
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
                        new_cost = cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
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
        else
            for n1=1:length(smallest_tour)
                if n1 == length(smallest_tour)
                    n3 = 1
                else
                    n3 = n1+1
                end
                cost1 = - T[smallest_tour[n1]+1, smallest_tour[n3]+1]
                for (i,tour) in enumerate(subtours)
                    for j=1:length(tour)                    
                        if j == length(tour)
                            n2 = j
                            n4 = 1
                        else
                            n2 = j
                            n4 = j+1
                        end
                        cost = cost1 - T[tour[n2]+1, tour[n4]+1]
                        new_cost = cost + T[smallest_tour[n1]+1, tour[n2]+1] + T[tour[n4]+1, smallest_tour[n3]+1]
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
        end
#         println("n1 = $best_n1")
#         println("n2 = $best_n2")
#         println("n3 = $best_n3")
#         println("n4 = $best_n4")
        if best_n3 == 1
            if best_n4 == 1 && best_n2 == length(subtours[best_tour])
                new_tour = vcat(smallest_tour, reverse(subtours[best_tour]))
            elseif best_n2 == 1 && best_n4 == length(subtours[best_tour])
                new_tour = vcat(smallest_tour, subtours[best_tour])
            else
                if best_n2 < best_n4
                    new_tour = vcat(reverse(smallest_tour), subtours[best_tour][best_n4:end], subtours[best_tour][1:best_n2])
                else
                    new_tour = vcat(smallest_tour, subtours[best_tour][best_n2:end], subtours[best_tour][1:best_n4])
                end
            end
        else
            if best_n4 == 1 && best_n2 == length(subtours[best_tour])
                new_tour = vcat(smallest_tour[1:best_n1], reverse(subtours[best_tour]), smallest_tour[best_n3:end])
            elseif best_n2 == 1 && best_n4 == length(subtours[best_tour])
                new_tour = vcat(smallest_tour[1:best_n1], subtours[best_tour], smallest_tour[best_n3:end])
            else
                if best_n2 < best_n4
                    new_tour = vcat(smallest_tour[1:best_n1], reverse(subtours[best_tour][1:best_n2]), reverse(subtours[best_tour][best_n4:end]), smallest_tour[best_n3:end])
                else
                    new_tour = vcat(smallest_tour[1:best_n1], subtours[best_tour][best_n2:end], subtours[best_tour][1:best_n4], smallest_tour[best_n3:end])
                end
            end
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

function eax_rand(p1::Vector{Int}, p2::Vector{Int}, T::Matrix{Float64})
    n = length(p1)
    AB_Cycles, EA = find_AB_cycles(p1, p2, n)
    select_effective_cycles(AB_Cycles)
    Eset = form_Eset_rand(AB_Cycles)
    solutions = generate_intermediate_solution(EA, Eset)
    subtours = separate_subtours(solutions, n)
    new_tour = modify_intermediate_solution(subtours, T)
    return new_tour
end

function generate_subtours_fast(tourA::Vector{Int}, Eset::Vector{abEdge})
    city = copy(tourA)
    push!(city, 0)
    n = length(city)
    pos = OffsetArray(zeros(Int, n), 0:n-1)
    # pos = zeros(Int, n)
    for i=1:n
        pos[city[i]] = i
    end

    dashed = Float64[]
    connectors = Vector{Vector{Int}}()
    for edge in Eset
        if edge.first
            if pos[edge.edge[1]] == n && pos[edge.edge[2]] == 1
                push!(dashed, (2*n+1)/2)
            elseif pos[edge.edge[2]] == n && pos[edge.edge[1]] == 1
                push!(dashed, (2*n+1)/2)
            else
                push!(dashed, (pos[edge.edge[1]]+pos[edge.edge[2]])/2)
            end
        else
            if pos[edge.edge[1]] < pos[edge.edge[2]]
                push!(connectors, [pos[edge.edge[1]], pos[edge.edge[2]]])
            else
                push!(connectors, [pos[edge.edge[2]], pos[edge.edge[1]]])
            end
        end
    end

    sort!(dashed)
    push!(dashed, dashed[1])

    segments = Vector{Vector{Int}}()
    for i=1:length(dashed)-1
        push!(segments, [Int(ceil(dashed[i])), Int(floor(dashed[i+1]))])
    end
    if segments[end][1] > n
        segments[length(segments)][1] = 1
    end
    mutual_segments = intersect(segments, connectors)
    subtours = Vector{Vector{Int}}()
    for segment in mutual_segments
        subtour = Int[]
        for i=segment[1]:segment[2]
            push!(subtour, city[i])
        end
        push!(subtours, subtour)
    end

    setdiff!(segments, mutual_segments)
    setdiff!(connectors, mutual_segments)
    sort!(segments, by=x->x[1])
    k = findfirst(x->x==reverse(segments[end]), connectors) 
    

    if !isnothing(k)
        subtour = Int[]
        deleteat!(connectors, k)
        b, a = segments[end]
        for i=0:n+a-b
            p = i+b
            if p <= n
                push!(subtour, city[p])
            else
                push!(subtour, city[p-n])
            end
        end
        push!(subtours, subtour)
        pop!(segments)
    end
    while !isempty(connectors)
        subtour = Int[]
        a, b = segments[1]
        for i = a:b
            push!(subtour,city[i])
        end
        deleteat!(segments, 1)
        while true
            i = findfirst(x->x[1]==b, connectors)
            if isnothing(i)
                j = findfirst(x->x[2]==b, connectors)
                c = connectors[j][1]
                deleteat!(connectors, j)
            else
                c = connectors[i][2]
                deleteat!(connectors, i)
            end
            if c == a 
                push!(subtours, subtour)
                break
            end
            k = findfirst(x->x[1]==c, segments)
            if isnothing(k)
                k = findfirst(x->x[2]==c, segments)
                d = segments[k][1]
                if segments[k][1] > segments[k][2]
                    for i = 0:n-d+c
                        p = c - i
                        if p > 0
                            push!(subtour, city[p])
                        else
                            push!(subtour, city[p+n])
                        end
                    end
                else
                    for i in reverse([l for l=segments[k][1]:segments[k][2]])
                        push!(subtour, city[i])
                    end
                end
            else
                d = segments[k][2]
                if segments[k][1] > segments[k][2]
                    for i = 0:n+d-c
                        p = c + i
                        if p <= n
                            push!(subtour, city[p])
                        else
                            push!(subtour, city[p-n])
                        end
                    end
                else
                    for i in [l for l=segments[k][1]:segments[k][2]]
                        push!(subtour, city[i])
                    end
                end
            end
            deleteat!(segments, k)
            b = d
        end
    end
    return subtours
end

function eax_1ab(p1::Vector{Int}, p2::Vector{Int}, T::Matrix{Float64})
    n = length(p1)
    AB_Cycles, EA = find_AB_cycles(p1, p2, n)
    select_effective_cycles(AB_Cycles)
    Eset = form_eset_1ab(AB_Cycles)
    subtours = generate_subtours_fast(p1, Eset)
    new_tour = modify_intermediate_solution(subtours, T)
    return new_tour
end

function have_intersects(a::Vector{Int}, b::Vector{Int})
    return any(x -> x in b, a)
end

function eax_block(p1::Vector{Int}, p2::Vector{Int}, T::Matrix{Float64}, nchr::Int)
    n = length(p1)
    AB_Cycles, EA = find_AB_cycles(p1, p2, n)
    select_effective_cycles(AB_Cycles)
    new_tours = Vector{Vector{Int}}()
    for k=1:min(nchr, length(AB_Cycles))
        Eset = AB_Cycles[k].edges
        subtours = generate_subtours_fast(p1, Eset)
        sort!(subtours, by=x->length(x), rev=true)
        for j = k+1:length(AB_Cycles)
            for i = 2:length(subtours)
                if have_intersects(subtours[i], AB_Cycles[j].nodes)
                    Eset = vcat(Eset, AB_Cycles[j].edges)
                    break
                end
            end
        end
        subtours = generate_subtours_fast(p1, Eset)
        new_tour = modify_intermediate_solution(subtours, T)
        push!(new_tours, new_tour)
    end
    return new_tours
end
