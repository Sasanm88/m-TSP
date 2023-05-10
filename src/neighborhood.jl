# function N1!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Int}, n_nodes::Int)   #Shift(0,1)
#     r1 = argmax([chrm.tours[i].cost for i=1:length(chrm.tours)])
#     routes = [i for i=1:length(chrm.tours)]
#     tour1 = chrm.tours[r1].sequence
#     cost1 = chrm.tours[r1].cost
#     for r2 in setdiff(routes, r1)    
#         tour2 = chrm.tours[r2].sequence
#         cost2 = chrm.tours[r2].cost
#         k1 = rand(1:length(tour1))
#         city1 = tour1[k1]
#         candidates = Int[] 
#         nt = length(tour2)
#         if nt == 0
#             candidates = [1]
#         elseif nt == 1
#             candidates = [1,2]
#         elseif nt == 2
#             candidates = [1,2,3]
#         else
#             if city1 in close_nodes[n_nodes+1,:] || tour2[1] in close_nodes[city1,:] 
#                 push!(candidates, 1)
#             end
#             for i=2:nt
#                 if tour2[i-1] in close_nodes[city1,:] || tour2[i] in close_nodes[city1,:]
#                     push!(candidates, i)
#                 end
#             end
#             if city1 in close_nodes[n_nodes+1,:] || tour2[nt] in close_nodes[city1,:] 
#                 push!(candidates, nt+1)
#             end
#         end

#         if length(candidates) > 0
#             k2 = candidates[rand(1:length(candidates))]

#             new_cost2 = Calculate_new_cost_add_one(tour2, cost2, city1, k2, TT, n_nodes)
#             new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
#             if new_cost2 - cost2 < 2 * (cost1 - new_cost1)

#                 if new_cost2 < cost1

#                     insert!(tour2, k2, city1)

#                     deleteat!(tour1, k1)
#                     chrm.tours[r1].cost = new_cost1
#                     chrm.tours[r2].cost = new_cost2
#                     chrm.genes = Int[]
#                     chrm.fitness = maximum([chrm.tours[i].cost for i=1:length(chrm.tours)])
#                     for tour in chrm.tours
#                         chrm.genes = vcat(chrm.genes, tour.sequence)
#                     end
#                     return 
#                 end
#             end
#         end
#     end
# end

function N1!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Shift(0,1)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    candidates = Int[]
    nt = length(tour2)
    if nt == 0
        candidates = [1]
    elseif nt == 1
        candidates = [1, 2]
    elseif nt == 2
        candidates = [1, 2, 3]
    else
        if close_nodes[n_nodes+1, city1] || close_nodes[city1, tour2[1]]
            push!(candidates, 1)
        end
        for i = 2:nt
            if close_nodes[city1, tour2[i-1]] || close_nodes[city1, tour2[i]]
                push!(candidates, i)
            end
        end
        if close_nodes[n_nodes+1, city1] || close_nodes[city1, tour2[nt]] 
            push!(candidates, nt + 1)
        end
    end

    if length(candidates) == 0
        return
    end
    k2 = candidates[rand(1:length(candidates))]

    new_cost2 = calculate_new_cost_add_one(tour2, cost2, city1, k2, TT, n_nodes)
    new_cost1 = calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
    #     if new_cost2 - cost2 > 2 * (cost1 - new_cost1)
    #         return 
    #     end
    if new_cost2 >= cost1
        return
    end

    insert!(tour2, k2, city1)

    deleteat!(tour1, k1)
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function N2!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Int}, n_nodes::Int)   #Swap(1,1)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    if length(tour2) == 0
        return
    end
    candidates = Int[]
    if length(tour2) == 1
        if close_nodes[n_nodes+1, city1]
            push!(candidates, 1)
        end
    else
        if close_nodes[n_nodes+1, city1] || close_nodes[city1, tour2[2]]
            push!(candidates, 1)
        end
        for i = 2:length(tour2)-1
            if close_nodes[city1, tour2[i-1]] || close_nodes[city1, tour2[i+1]]
                push!(candidates, i)
            end
        end
        if close_nodes[n_nodes+1, city1] || close_nodes[city1, tour2[end]]
            push!(candidates, length(tour2))
        end
    end
    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    k2 = candidates[rand(1:length(candidates))]
    city2 = tour2[k2]

    new_cost1, new_cost2 = calculate_new_cost_swap_one(tour1, cost1, city1, k1, tour2, cost2, city2, k2, TT, n_nodes)

    if new_cost1 >= cost1 || new_cost2 >= cost1
        return
    end
    tour1[k1] = city2
    tour2[k2] = city1
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function N3!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Shift(0,2)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    if length(tour1) < 2
        return
    end
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    candidates = Int[]
    if length(tour2) == 0
        candidates = [1]
    elseif length(tour2) == 1
        if close_nodes[n_nodes+1, city1] || close_nodes[city2, tour2[1]]
            push!(candidates, 1)
        end
        if close_nodes[n_nodes+1, city2] || close_nodes[city1, tour2[1]]
            push!(candidates, 2)
        end
    elseif length(tour2) == 2
        if close_nodes[n_nodes+1, city1] || close_nodes[city2, tour2[1]]
            push!(candidates, 1)
        end
        if close_nodes[city1, tour2[1]] || close_nodes[city2, tour2[2]]
            push!(candidates, 2)
        end
        if close_nodes[n_nodes+1, city2] || close_nodes[city1, tour2[2]]
            push!(candidates, 3)
        end
    else
        if close_nodes[n_nodes+1, city1] || close_nodes[city2, tour2[1]]
            push!(candidates, 1)
        end
        for i = 2:length(tour2)-1
            if close_nodes[city1, tour2[i-1]] || close_nodes[city2, tour2[i]]
                push!(candidates, i + 1)
            end
        end
        if close_nodes[n_nodes+1, city2] || close_nodes[city1, tour2[length(tour2)]]
            push!(candidates, length(tour2) + 1)
        end
    end
    if length(candidates) == 0
        return
    end
    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:length(tour2)+1)

    new_cost2, straight = calculate_new_cost_add_two(tour2, cost2, city1, city2, k2, TT, n_nodes)

    if new_cost2 >= cost1
        return
    end
    if straight
        insert!(tour2, k2, city1)
        insert!(tour2, k2 + 1, city2)
    else
        insert!(tour2, k2, city2)
        insert!(tour2, k2 + 1, city1)
    end
    new_cost1 = calculate_new_cost_remove_two(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1 + 1])
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function N4!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Swap(2,2)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    if length(tour1) < 2
        return
    end
    k1 = rand(1:length(tour1)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    if length(tour2) < 2
        return
    end
    nt2 = length(tour2)
    candidates = Int[]
    if nt2 == 2
        if close_nodes[n_nodes+1, city11] || close_nodes[n_nodes+1, city12]
            push!(candidates, 1)
        end
    elseif nt2 == 3
        if close_nodes[n_nodes+1, city11] || close_nodes[city12, tour2[3]]
            push!(candidates, 1)
        end
        if close_nodes[tour2[nt2-2], city11] || close_nodes[n_nodes+1, city12]
            push!(candidates, nt2 - 1)
        end
    else
        if close_nodes[n_nodes+1, city11] || close_nodes[city12, tour2[3]]
            push!(candidates, 1)
        end
        for i = 2:nt2-2
            if close_nodes[city11, tour2[i-1]] || close_nodes[city12, tour2[i+2]]
                push!(candidates, i)
            end
        end
        if close_nodes[tour2[nt2-2], city11] || close_nodes[n_nodes+1, city12]
            push!(candidates, nt2 - 1)
        end
    end
    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_swap_two_updated(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return
    end
    if straight2
        tour1[k1] = city21
        tour1[k1+1] = city22
    else
        tour1[k1] = city22
        tour1[k1+1] = city21
    end
    if straight1
        tour2[k2] = city11
        tour2[k2+1] = city12
    else
        tour2[k2] = city12
        tour2[k2+1] = city11
    end
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function N5!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Shift(0,3)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    if length(tour1) < 3
        return
    end
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-2)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    city3 = tour1[k1+2]
    candidates = Int[]
    if length(tour2) == 0
        candidates = [1]
    elseif length(tour2) == 1
        if close_nodes[n_nodes+1, city1] || close_nodes[city3, tour2[1]]
            push!(candidates, 1)
        end
        if close_nodes[n_nodes+1, city3] || close_nodes[city1, tour2[1]]
            push!(candidates, 2)
        end
    elseif length(tour2) == 3
        if close_nodes[n_nodes+1, city1] || close_nodes[city3, tour2[1]]
            push!(candidates, 1)
        end
        if close_nodes[city1, tour2[1]] || close_nodes[city3, tour2[2]]
            push!(candidates, 2)
        end
        if close_nodes[n_nodes+1, city3] || close_nodes[city1, tour2[2]]
            push!(candidates, 3)
        end
    else
        if close_nodes[n_nodes+1, city1] || close_nodes[city3, tour2[1]]
            push!(candidates, 1)
        end
        for i = 2:length(tour2)-1
            if close_nodes[city1, tour2[i-1]] || close_nodes[city3, tour2[i]]
                push!(candidates, i + 1)
            end
        end
        if close_nodes[n_nodes+1, city3] || close_nodes[city1, tour2[length(tour2)]]
            push!(candidates, length(tour2) + 1)
        end
    end
    if length(candidates) == 0
        return
    end
    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:length(tour2)+1)

    new_cost2, straight = calculate_new_cost_add_three(tour2, cost2, city1, city2, city3, k2, TT, n_nodes)

    if new_cost2 >= cost1
        return
    end
    if straight
        insert!(tour2, k2, city1)
        insert!(tour2, k2 + 1, city2)
        insert!(tour2, k2 + 2, city3)
    else
        insert!(tour2, k2, city3)
        insert!(tour2, k2 + 1, city2)
        insert!(tour2, k2 + 2, city1)
    end

    new_cost1 = calculate_new_cost_remove_three(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1 + 1, k1 + 2])
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end


function N6!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Swap(3,3)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 3
        return
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    candidates = Int[]
    if nt2 == 3
        if close_nodes[n_nodes+1, city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, 1)
        end
    elseif nt2 == 4  #Continue from here
        if close_nodes[n_nodes+1, city11] || close_nodes[city13, tour2[4]]
            push!(candidates, 1)
        end
        if close_nodes[tour2[nt2-3], city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, 2)
        end
    else
        if close_nodes[n_nodes+1, city11] || close_nodes[city13, tour2[4]]
            push!(candidates, 1)
        end
        for i = 2:nt2-3
            if close_nodes[city11, tour2[i-1]] || close_nodes[city13, tour2[i+3]]
                push!(candidates, i)
            end
        end
        if close_nodes[tour2[nt2-3], city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, nt2 - 2)
        end
    end
    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    k2 = candidates[rand(1:length(candidates))]
    #     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]

    new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_swap_three_updated(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return
    end
    if straight2
        tour1[k1] = city21
        tour1[k1+1] = city22
        tour1[k1+2] = city23
    else
        tour1[k1] = city23
        tour1[k1+1] = city22
        tour1[k1+2] = city21
    end
    if straight1
        tour2[k2] = city11
        tour2[k2+1] = city12
        tour2[k2+2] = city13
    else
        tour2[k2] = city13
        tour2[k2+1] = city12
        tour2[k2+2] = city11
    end
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function N7!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int)   #Swap(3,2)
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    tour1 = chrm.tours[r1].sequence
    tour2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 2
        return
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    candidates = Int[]
    if nt2 == 2
        if close_nodes[n_nodes+1, city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, 1)
        end
    elseif nt2 == 3
        if close_nodes[n_nodes+1, city11] || close_nodes[city13, tour2[3]]
            push!(candidates, 1)
        end
        if close_nodes[tour2[1], city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, 2)
        end
    else
        if close_nodes[n_nodes+1, city11] || close_nodes[city13, tour2[3]]
            push!(candidates, 1)
        end
        for i = 2:nt2-2
            if close_nodes[city11, tour2[i-1]] || close_nodes[city13, tour2[i+2]]
                push!(candidates, i)
            end
        end
        if close_nodes[tour2[nt2-2], city11] || close_nodes[n_nodes+1, city13]
            push!(candidates, nt2 - 1)
        end
    end
    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    k2 = candidates[rand(1:length(candidates))]
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_swap_three_with_two_updated(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return
    end
    if straight2
        tour1[k1] = city21
        tour1[k1+1] = city22
        deleteat!(tour1, k1 + 2)
    else
        tour1[k1] = city22
        tour1[k1+1] = city21
        deleteat!(tour1, k1 + 2)
    end
    if straight1
        tour2[k2] = city11
        tour2[k2+1] = city12
        insert!(tour2, k2 + 2, city13)
    else
        tour2[k2] = city13
        tour2[k2+1] = city12
        insert!(tour2, k2 + 2, city11)
    end
    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end


function N_cross!(chrm::Chromosome, T::Matrix{Float64}, close_nodes::Matrix{Int}, n_nodes::Int)   #Cross Exchange
    r1 = argmax([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    routes = [i for i in 1:length(chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(chrm.tours)-1)]
    t1 = chrm.tours[r1].sequence
    t2 = chrm.tours[r2].sequence
    cost1 = chrm.tours[r1].cost
    cost2 = chrm.tours[r2].cost
    if length(t1) < 6 || length(t2) < 6
        return
    end
    tau = min(8, Int(round(min(length(t1), length(t2)) / 2)))
    k11 = rand(1:length(t1)-tau)
    l1 = rand(1:tau)
    k12 = k11 + l1

    candidates = Int[]
    for i in 1:length(t2)-tau
        if k11 == 1
            if close_nodes[n_nodes+1, t2[i]]
                push!(candidates, 1)
            end
        else
            if close_nodes[t1[k11-1], t2[i]]
                push!(candidates, i)
            end
        end
    end

    if length(candidates) == 0
        return
    end
    candidates = collect(Set(candidates))
    k21 = candidates[rand(1:length(candidates))]
    l2 = rand(1:tau)
    k22 = k21 + l2

    new_cost1, new_cost2, straight1, straight2 = calculate_new_cost_cross(t1, cost1, t2, cost2, k11, k12, k21, k22, T, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return
    end

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
    deleteat!(t1, [i for i = k11:k12])
    for i in 1:k22-k21+1
        insert!(t1, i + k11 - 1, alpha2[i])
    end

    deleteat!(t2, [i for i = k21:k22])
    for i in 1:k12-k11+1
        insert!(t2, i + k21 - 1, alpha1[i])
    end

    chrm.tours[r1].cost = new_cost1
    chrm.tours[r2].cost = new_cost2
    chrm.fitness = maximum([chrm.tours[i].cost for i in 1:length(chrm.tours)])
    index = 0
    for tour in chrm.tours
        chrm.genes[index+1:index+length(tour.sequence)] = tour.sequence
        index += length(tour.sequence)
    end
    return
end

function Improve_chromosome!(chrm::Chromosome, TT::Matrix{Float64}, close_nodes::Matrix{Bool}, n_nodes::Int, roullet::Vector{Int}, old_best::Float64, iterations::Int)

    for i in 1:iterations
        r = sample(1:length(roullet), weights(roullet))

        f1::Float64 = chrm.fitness

        if r == 1
            Ni1!(chrm, TT, close_nodes, n_nodes)
        elseif r == 2
            Ni2!(chrm, TT, close_nodes, n_nodes)
        elseif r == 3
            Ni3!(chrm, TT, close_nodes, n_nodes)
        elseif r == 4
            Ni4!(chrm, TT, close_nodes, n_nodes)
        elseif r == 5
            Ni5!(chrm, TT, close_nodes, n_nodes)
        elseif r == 6
            N1!(chrm, TT, close_nodes, n_nodes)
        else
            error("This should not happen...")
        end

        if chrm.fitness < f1
            roullet[r] += 1
        end
    end

end

function Improve_Population!(P::Vector{Chromosome}, TT::Matrix{Float64}, close_nodes::Matrix{Int}, n_nodes::Int)
    #     Search_methods = [N4, N2, N3, Ni6, Ni7, N3r, N4sr, N4rs, N4rr, N5, N5r, N6, N6sr, N6rs, N6rr, N7, N7rs, N7sr, N7rr]
    #     Search_methods = [N1, N2, N3, N4, N5, N6, N7, Ni1, Ni2, Ni3, Ni4, Ni5, Ni6, Ni7, N_cross]
    #     seq = sort(sample(1:length(P), 10, replace = false))[1:5]
    #     seq = [1]
    #     for j in seq
    #         chrm = P[j]
    #     Search_methods = Function[N1!, Ni1!, Ni2!, Ni3!, Ni4!, Ni5!, Ni6!, Ni7!]
    #     for chrm in P
    #         for i=1:1000
    #             r = rand(1:length(Search_methods))
    #             search = Search_methods[r]
    #             f = chrm.fitness
    #             chrm = search(chrm, TT, close_nodes, n_nodes)
    # # 
    #         end
    #     end
    #     chrm = P[1]
    #     P[1] = N_cross_extensive(chrm, T, n_nodes, 5)
end

