using StatsBase


function Crossover_OX1(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #order crossover
    child = zeros(Int64, n_nodes)
    idx1 = rand(2:n_nodes-1)
    idx2 = rand(2:n_nodes-1)
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    copy_from_P1 = parent1[idx1:idx2] 
    child[idx1:idx2] = copy_from_P1
    index = 1
    for gene in parent2
        if index == idx1
            if idx2 == n_nodes
                break
            else
                index = idx2+1
            end
        end
        if !(gene in copy_from_P1)
            child[index] = gene
            index += 1
        end
    end
    return child
end

function Crossover_POS(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)  #position based crossover
    child = zeros(Int64, n_nodes)
    num_pos = rand(1:n_nodes-1)
    selected_pos = sample(1:n_nodes,Weights(ones(n_nodes)),num_pos, replace=false)
    selected_p1 = parent1[selected_pos]
    child[selected_pos] = selected_p1

    for i in parent2
        if !(i in selected_p1)
            child[findfirst(x->x==0, child)] = i
        end
    end
    return child
end


function Crossover_CX(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)  #Cyclic crossover 
    child = zeros(Int64, n_nodes)
    cycles = []
    labels = zeros(n_nodes)
    abs_p2 = abs.(parent2)
    while 0 in labels
        cycle = []
        idx1 = findfirst(x->x == 0, labels)
        labels[idx1] = 1
        node1 = parent1[idx1]
        push!(cycle, idx1)
        while true
            idx = findfirst(x->x == node1, parent2)
            if idx == idx1
                break
            else
                push!(cycle,idx)
                node1 = parent1[idx]
                labels[idx] = 1
            end
        end
        push!(cycles,cycle)
    end
    for (i,cycle) in enumerate(cycles)
        for c in cycle
            if i%2==0
                child[c]=parent1[c]
            else
                child[c] = parent2[c]
            end
        end
    end
    return child
end


function Crossover_OX2(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #Order-based crossover
    child = zeros(Int64, n_nodes)
    num_pos = rand(1:n_nodes)
    selected_pos2 = sample(1:n_nodes,Weights(ones(n_nodes)),num_pos, replace=false)
    selected_p2 = parent2[selected_pos2]
    
    selected_pos1 = findall(x->x in selected_p2, parent1)
    unselected_pos = setdiff(1:n_nodes,selected_pos1)
    child[unselected_pos] = parent1[unselected_pos]
    child[selected_pos1] = parent2[sort(selected_pos2)]
    return child
end


function find_mapping(Dict, parent_selected_part, i)
    while true 
        j = Dict[i]
        if !(j in parent_selected_part)
            return j
        else
            i=j
        end
    end
end

function Crossover_PMX(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #Partially mapped crossover
    child1 = zeros(Int64, n_nodes)
#     child2 = zeros(Int64, n_nodes)
    idx1 = rand(2:n_nodes-1)
    idx2 = rand(2:n_nodes-1)
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end

    copy_from_P1 = parent1[idx1:idx2] 
    copy_from_P2 = parent2[idx1:idx2]
    child1[idx1:idx2] = copy_from_P2
    child1[1:idx1-1] = parent1[1:idx1-1]
    child1[idx2+1:n_nodes] = parent1[idx2+1:n_nodes]
#     child2[idx1:idx2] = copy_from_P1
#     child2[1:idx1-1] = parent2[1:idx1-1]
#     child2[idx2+1:n_nodes] = parent2[idx2+1:n_nodes]
#     change1 = Dict{Int, Int}()
    change2 = Dict{Int, Int}()
    for i=1:idx2-idx1+1
#         change1[copy_from_P1[i]] = copy_from_P2[i]
        change2[copy_from_P2[i]] = copy_from_P1[i]
    end
    for i=1:n_nodes
        if i < idx1 || i>idx2
            if child1[i] in copy_from_P2
                child1[i] = find_mapping(change2, copy_from_P2, child1[i])
            end
#             if child2[i] in copy_from_P1
#                 child2[i] = find_mapping(change1, copy_from_P1, child2[i])
#             end
        end
    end
    return child1 #, child2
end


function Crossover_HX(TT::Matrix{Float64}, parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #heuristic crossover 
    
    remaining_cities = copy(parent1)
    current_city = rand(1:n_nodes)
    child = [current_city]
    deleteat!(remaining_cities, findfirst(x->x==current_city, remaining_cities))
    while length(remaining_cities)>0
        pos1 = findfirst(x->x==current_city, parent1)
        pos2 = findfirst(x->x==current_city, parent2)
        next_city = n_nodes+2
        min_edge = 100000

        if pos1 == 1
            if !(parent1[pos1+1] in child)
                if TT[current_city+1,parent1[pos1+1]+1] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city+1,parent1[pos1+1]+1]
                end
            end
        elseif pos1==n_nodes
            if !(parent1[pos1-1] in child)
                if TT[current_city+1,parent1[pos1-1]+1] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city+1,parent1[pos1-1]+1]
                end
            end
        else
            if !(parent1[pos1+1] in child)
                if TT[current_city+1,parent1[pos1+1]+1] < min_edge
                    next_city = parent1[pos1+1]
                    min_edge = TT[current_city+1,parent1[pos1+1]+1]
                end
            end
            if !(parent1[pos1-1] in child)
                if TT[current_city+1,parent1[pos1-1]+1] < min_edge
                    next_city = parent1[pos1-1]
                    min_edge = TT[current_city+1,parent1[pos1-1]+1]
                end
            end
        end

        if pos2 == 1
            if !(parent2[pos2+1] in child)
                if TT[current_city+1, parent2[pos2+1]+1] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city+1,parent2[pos2+1]+1]
                end
            end
        elseif pos2==n_nodes
            if !(parent2[pos2-1] in child)
                if TT[current_city+1, parent2[pos2-1]+1] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city+1,parent2[pos2-1]+1]
                end
            end
        else
            if !(parent2[pos2+1] in child)
                if TT[current_city+1, parent2[pos2+1]+1] < min_edge
                    next_city = parent2[pos2+1]
                    min_edge = TT[current_city+1,parent2[pos2+1]+1]
                end
            end
            if !(parent2[pos2-1] in child)
                if TT[current_city+1, parent2[pos2-1]+1] < min_edge
                    next_city = parent2[pos2-1]
                    min_edge = TT[current_city+1,parent2[pos2-1]+1]
                end
            end
        end

        if next_city == n_nodes+2
            next_city = remaining_cities[rand(1:length(remaining_cities))]
        end

        current_city = next_city
        push!(child, next_city)
        deleteat!(remaining_cities, findfirst(x->x==current_city, remaining_cities))
    end
    return child
end


function Remove_one_city(tours::Vector{Tour}, city::Int, T::Matrix{Float64}, n_nodes::Int)
    for (t,tour) in enumerate(tours)
        nt = length(tour.Sequence)
        k = findfirst(x->x==city, tour.Sequence)
        if !isnothing(k)
            if k == 1
                if nt == 1
                    deleteat!(tours, t)
                else
                    tour.cost = tour.cost - T[1, city+1] - T[city+1, tour.Sequence[2]+1] + T[1, tour.Sequence[2]+1]
                    deleteat!(tour.Sequence, 1)
                end
            elseif k == nt
                tour.cost = tour.cost - T[city+1, n_nodes+2] - T[tour.Sequence[nt-1]+1, city+1] + T[tour.Sequence[nt-1]+1, n_nodes+2]
                deleteat!(tour.Sequence, nt)
            else
                tour.cost = tour.cost - T[tour.Sequence[k-1]+1, city+1] - T[city+1, tour.Sequence[k+1]+1] + T[tour.Sequence[k-1]+1, tour.Sequence[k+1]+1]
                deleteat!(tour.Sequence, k)
            end
        end
    end
end

function put_city_in_tour(c::Vector{Tour}, city::Int, T::Matrix{Float64}, n_nodes::Int)
    least_increase = Inf
    best_tour = 0
    best_position = 0
    for i=2:length(c)
        tour = c[i].Sequence
        nt = length(tour)
        if nt==0
            increase = T[1, city+1] + T[city+1, n_nodes+2]
            best_tour = i
            best_position = 1
            break
        end
        increase = T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1]
        if increase < least_increase
            least_increase = increase
            best_tour = i
            best_position = 1
        end
        for j = 2:nt
            increase = T[tour[j-1]+1, city+1] + T[city+1, tour[j]+1] - T[tour[j-1]+1, tour[j]+1]
            if increase < least_increase
                least_increase = increase
                best_tour = i
                best_position = j
            end
        end
        increase = T[tour[nt]+1, city+1] + T[city+1, n_nodes+2] - T[tour[nt]+1, n_nodes+2]
        if increase < least_increase
            least_increase = increase
            best_tour = i
            best_position = nt+1
        end
    end
    insert!(c[best_tour].Sequence, best_position, city)
    c[best_tour].cost += least_increase
    if c[best_tour].cost > c[1].cost
        temp = deepcopy(c[1])
        c[1] = c[best_tour]
        c[best_tour] = temp
    end
end          

function new_crossover(parent1::Chromosome, parent2::Chromosome, T::Matrix{Float64}, n_nodes::Int64)
    P1 = deepcopy(parent1)
    P2 = deepcopy(parent2)
    if length(P1.tours) == 0
        println(P1)
    end
    c = Tour[]
    max_tour_length = 0.0
#     m = min(length(P1.tours), length(P2.tours)) - 1
    m = length(P1.tours)
    for i=1:m
        if length(P1.tours) + length(P2.tours) == 0
            break
        end
        if length(P2.tours) == 0
            pr = 1
        elseif length(P1.tours) == 0
            pr = 2
        else
            if rand() < 0.5
                pr = 1
            else
                pr = 2
            end
        end
        
        if pr == 1
            if length(P1.tours) + length(P2.tours) == 0
                println("It Happened")
            end
            r = argmin([tour.cost for tour in P1.tours])
#             println("tour ", r, " chosen from parent 1")
            if P1.tours[r].cost > max_tour_length
                max_tour_length = P1.tours[r].cost
                pushfirst!(c, P1.tours[r])
            else
                push!(c, P1.tours[r])
            end
            for city in P1.tours[r].Sequence
                Remove_one_city(P2.tours, city, T, n_nodes)
            end
            deleteat!(P1.tours, r)
        else
            r = argmin([tour.cost for tour in P2.tours])
#             println("tour ", r, " chosen from parent 2")
            if P2.tours[r].cost > max_tour_length
                max_tour_length = P2.tours[r].cost
                pushfirst!(c, P2.tours[r])
            else
                push!(c, P2.tours[r])
            end
            
            for city in P2.tours[r].Sequence
                Remove_one_city(P1.tours, city, T, n_nodes)
            end
            deleteat!(P2.tours, r)
        end
    end
    Remaining = Int[]
    for tour in P1.tours
        for city in tour.Sequence
            push!(Remaining, city)
        end
    end
    for tour in P2.tours
        for city in tour.Sequence
            push!(Remaining, city)
        end
    end
    R = collect(Set(Remaining))
    for city in R 
        put_city_in_tour(c, city, T, n_nodes)
    end
    child = Int[]
    for tour in c
        for city in tour.Sequence
            push!(child, city)
        end
    end
    return child
end