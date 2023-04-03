
function N1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,1)

    r1, r2 = sample(1:length(Chrm.tours), 2, replace = false)
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) == 1
        return Chrm
    end
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    Candidates = Int[] 
    nt = length(tour2)
    if nt == 1
        Candidates = [1,2]
    elseif nt == 2
        Candidates = [1,2,3]
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city1,:] 
            push!(Candidates, 1)
        end
        for i=2:nt
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i] in Close_nodes[city1,:]
                push!(Candidates, i)
            end
        end
        if city1 in Close_nodes[n_nodes+1,:] || tour2[nt] in Close_nodes[city1,:] 
            push!(Candidates, nt+1)
        end
    end
    # Candidates = collect(setdiff(Set(Candidates),Set([k1])))
    
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)+1)
    
    new_cost2 = Calculate_new_cost_add_one(tour2, cost2, city1, k2, TT, n_nodes)
    new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
    
    if new_cost2 + new_cost1 >= cost1 + cost2
        return Chrm
    end
    insert!(tour2, k2, city1)
    
    deleteat!(tour1, k1)
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = Chrm.fitness + new_cost1 + new_cost2 - cost1 - cost2
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end


function N2(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(1,1)
    r1, r2 = sample(1:length(Chrm.tours), 2, replace = false)
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    if length(tour2) == 0
        return Chrm
    end
    Candidates = Int[] 
    if length(tour2) == 1
        if city1 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[2] in Close_nodes[city1,:] 
            push!(Candidates, 1)
        end
        for i=2:length(tour2)-1
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i+1] in Close_nodes[city1,:]
                push!(Candidates, i)
            end
        end
        if city1 in Close_nodes[n_nodes+1,:] || tour2[length(tour2)] in Close_nodes[city1,:] 
            push!(Candidates, length(tour2))
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
    city2 = tour2[k2]

    new_cost1, new_cost2 = Calculate_new_cost_swap_one(tour1, cost1, city1, k1, tour2, cost2, city2, k2, TT, n_nodes)

    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city2
    tour2[k2] = city1
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N3(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,2)
    r1, r2 = sample(1:length(Chrm.tours), 2, replace = false)
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) < 2
        return Chrm
    end
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    Candidates = Int[] 
    if length(tour2) == 0
        Candidates = [1]
    elseif length(tour2) == 1
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city1,:] 
            push!(Candidates, 2)
        end
    elseif length(tour2) == 2
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        if tour2[1] in Close_nodes[city1,:] || tour2[2] in Close_nodes[city2,:]
            push!(Candidates, 2)
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[2] in Close_nodes[city1,:] 
            push!(Candidates, 3)
        end
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        for i=2:length(tour2)-1
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i] in Close_nodes[city2,:]
                push!(Candidates, i+1)
            end
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[length(tour2)] in Close_nodes[city1,:] 
            push!(Candidates, length(tour2)+1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)+1)

    new_cost2 = Calculate_new_cost_add_two(tour2, cost2, city1, city2, k2, TT, n_nodes)
    
    if new_cost2 >= cost1
        return Chrm
    end

    insert!(tour2, k2, city1)
    insert!(tour2, k2+1, city2)
    
    new_cost1 = Calculate_new_cost_remove_two(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1+1])
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    return Chrm
end

function N3r(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) < 2
        return Chrm
    end
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    Candidates = Int[] 
    if length(tour2) == 0
        Candidates = [1]
    elseif length(tour2) == 1
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city1,:] 
            push!(Candidates, 2)
        end
    elseif length(tour2) == 2
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        if tour2[1] in Close_nodes[city1,:] || tour2[2] in Close_nodes[city2,:]
            push!(Candidates, 2)
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[2] in Close_nodes[city1,:] 
            push!(Candidates, 3)
        end
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        for i=2:length(tour2)-1
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i] in Close_nodes[city2,:]
                push!(Candidates, i+1)
            end
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour2[length(tour2)] in Close_nodes[city1,:] 
            push!(Candidates, length(tour2)+1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)+1)

    new_cost2 = Calculate_new_cost_add_two(tour2, cost2, city2, city1, k2, TT, n_nodes)
    
    if new_cost2 >= cost1
        return Chrm
    end

    insert!(tour2, k2, city2)
    insert!(tour2, k2+1, city1)
    
    new_cost1 = Calculate_new_cost_remove_two(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1+1])
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    return Chrm
end


function N4(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 2 
        return Chrm
    end
    k1 = rand(1:length(tour1)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    if length(tour2) < 2 
        return Chrm
    end
    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2
        if city11 in Close_nodes[n_nodes+1,:] || city12 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+2] in Close_nodes[city12,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    tour2[k2] = city11
    tour2[k2+1] = city12
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end


function N4rs(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 2 
        return Chrm
    end
    k1 = rand(1:length(tour1)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    if length(tour2) < 2 
        return Chrm
    end
    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2
        if city11 in Close_nodes[n_nodes+1,:] || city12 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+2] in Close_nodes[city12,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two_reverse_straight(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    tour2[k2] = city12
    tour2[k2+1] = city11
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N4sr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour2) < 2 
        return Chrm
    end
    k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    if length(tour1) < 2 
        return Chrm
    end
    nt1 = length(tour1)
    Candidates = Int[] 
    if nt1 == 2
        if city21 in Close_nodes[n_nodes+1,:] || city22 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt1 == 3
        if city21 in Close_nodes[n_nodes+1,:] || tour1[3] in Close_nodes[city22,:] 
            push!(Candidates, 1)
        end
        if city21 in Close_nodes[tour1[nt1-2],:] || city22 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt1-1)
        end
    else
        if city21 in Close_nodes[n_nodes+1,:] || tour1[3] in Close_nodes[city22,:] 
            push!(Candidates, 1)
        end
        for i=2:nt1-2
            if tour1[i-1] in Close_nodes[city21,:] || tour1[i+2] in Close_nodes[city22,:]
                push!(Candidates, i)
            end
        end
        if city21 in Close_nodes[tour1[nt1-2],:] || city22 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt1-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k1 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two_straight_reverse(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city22
    tour1[k1+1] = city21
    tour2[k2] = city11
    tour2[k2+1] = city12
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N4rr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 2 
        return Chrm
    end
    k1 = rand(1:length(tour1)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    if length(tour2) < 2 
        return Chrm
    end
    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2
        if city11 in Close_nodes[n_nodes+1,:] || city12 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city12,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+2] in Close_nodes[city12,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city12 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two_reverse_reverse(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city22
    tour1[k1+1] = city21
    tour2[k2] = city12
    tour2[k2+1] = city11
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N5(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) < 3
        return Chrm
    end
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-2)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    city3 = tour1[k1+2]
    Candidates = Int[] 
    if length(tour2) == 0
        Candidates = [1]
    elseif length(tour2) == 1
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city1,:] 
            push!(Candidates, 2)
        end
    elseif length(tour2) == 3
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        if tour2[1] in Close_nodes[city1,:] || tour2[2] in Close_nodes[city3,:]
            push!(Candidates, 2)
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[2] in Close_nodes[city1,:] 
            push!(Candidates, 3)
        end
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        for i=2:length(tour2)-1
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i] in Close_nodes[city3,:]
                push!(Candidates, i+1)
            end
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[length(tour2)] in Close_nodes[city1,:] 
            push!(Candidates, length(tour2)+1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)+1)

    new_cost2 = Calculate_new_cost_add_three(tour2, cost2, city1, city2, city3, k2, TT, n_nodes)
    
    if new_cost2 >= cost1
        return Chrm
    end

    insert!(tour2, k2, city1)
    insert!(tour2, k2+1, city2)
    insert!(tour2, k2+2, city3)
    
    new_cost1 = Calculate_new_cost_remove_three(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1+1, k1+2])
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    return Chrm
end

function N5r(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) < 3
        return Chrm
    end
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    k1 = rand(1:length(tour1)-2)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    city3 = tour1[k1+2]
    Candidates = Int[] 
    if length(tour2) == 0
        Candidates = [1]
    elseif length(tour2) == 1   #Change this for reverse
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city1,:] 
            push!(Candidates, 2)
        end
    elseif length(tour2) == 3
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        if tour2[1] in Close_nodes[city1,:] || tour2[2] in Close_nodes[city3,:]
            push!(Candidates, 2)
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[2] in Close_nodes[city1,:] 
            push!(Candidates, 3)
        end
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour2[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        for i=2:length(tour2)-1
            if tour2[i-1] in Close_nodes[city1,:] || tour2[i] in Close_nodes[city3,:]
                push!(Candidates, i+1)
            end
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour2[length(tour2)] in Close_nodes[city1,:] 
            push!(Candidates, length(tour2)+1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)+1)

    new_cost2 = Calculate_new_cost_add_three(tour2, cost2, city3, city2, city1, k2, TT, n_nodes)
    
    if new_cost2 >= cost1
        return Chrm
    end

    insert!(tour2, k2, city3)
    insert!(tour2, k2+1, city2)
    insert!(tour2, k2+2, city1)
    
    new_cost1 = Calculate_new_cost_remove_three(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, [k1, k1+1, k1+2])
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    return Chrm
end

function N6(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 3
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 3  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  #Continue from here
        if city11 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[nt2-3],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-3
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+3] in Close_nodes[city13,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-3],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-2)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    tour1[k1+2] = city23
    tour2[k2] = city11
    tour2[k2+1] = city12
    tour2[k2+2] = city13
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N6rs(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]
    if length(tour2) < 3 
        return Chrm
    end
    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 3
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4
        if city11 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[nt2-3],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-2)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-3
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+3] in Close_nodes[city13,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-3],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-2)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]

    new_cost1, new_cost2 = Calculate_new_cost_swap_three_reverse_straight(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    tour1[k1+2] = city23
    tour2[k2] = city13
    tour2[k2+1] = city12
    tour2[k2+2] = city11
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N6sr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour2) < 3 
        return Chrm
    end
    k2 = rand(1:length(tour2)-2)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]
    if length(tour1) < 3 
        return Chrm
    end
    nt1 = length(tour1)
    Candidates = Int[] 
    if nt1 == 3
        if city21 in Close_nodes[n_nodes+1,:] || city23 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt1 == 4
        if city21 in Close_nodes[n_nodes+1,:] || tour1[4] in Close_nodes[city23,:] 
            push!(Candidates, 1)
        end
        if city21 in Close_nodes[tour1[nt1-3],:] || city23 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt1-2)
        end
    else
        if city21 in Close_nodes[n_nodes+1,:] || tour1[4] in Close_nodes[city23,:] 
            push!(Candidates, 1)
        end
        for i=2:nt1-3
            if tour1[i-1] in Close_nodes[city21,:] || tour1[i+3] in Close_nodes[city23,:]
                push!(Candidates, i)
            end
        end
        if city21 in Close_nodes[tour1[nt1-3],:] || city23 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt1-2)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k1 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    new_cost1, new_cost2 = Calculate_new_cost_swap_three_straight_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city23
    tour1[k1+1] = city22
    tour1[k1+2] = city21
    tour2[k2] = city11
    tour2[k2+1] = city12
    tour2[k2+2] = city13
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N6rr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 3
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 3  
        if city13 in Close_nodes[n_nodes+1,:] || city11 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  #Continue from here
        if city13 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        if city13 in Close_nodes[tour2[nt2-3],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city13 in Close_nodes[n_nodes+1,:] || tour2[4] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-3
            if tour2[i-1] in Close_nodes[city13,:] || tour2[i+3] in Close_nodes[city11,:]
                push!(Candidates, i)
            end
        end
        if city13 in Close_nodes[tour2[nt2-3],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-2)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_reverse_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city23
    tour1[k1+1] = city22
    tour1[k1+2] = city21
    tour2[k2] = city13
    tour2[k2+1] = city12
    tour2[k2+2] = city11
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N7(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 2
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3  
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[1],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+2] in Close_nodes[city13,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    deleteat!(tour1, k1+2)
    tour2[k2] = city11
    tour2[k2+1] = city12
    insert!(tour2, k2+2, city13)

    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N7rs(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 2
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2  
        if city13 in Close_nodes[n_nodes+1,:] || city11 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3  
        if city13 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        if city13 in Close_nodes[tour2[1],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city13 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city13,:] || tour2[i+2] in Close_nodes[city11,:]
                push!(Candidates, i)
            end
        end
        if city13 in Close_nodes[tour2[nt2-2],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_reverse_straight(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city21
    tour1[k1+1] = city22
    deleteat!(tour1, k1+2)
    tour2[k2] = city13
    tour2[k2+1] = city12
    insert!(tour2, k2+2, city11)

    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N7sr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 2
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3  
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        if city11 in Close_nodes[tour2[1],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city11 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city13,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city11,:] || tour2[i+2] in Close_nodes[city13,:]
                push!(Candidates, i)
            end
        end
        if city11 in Close_nodes[tour2[nt2-2],:] || city13 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_straight_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city22
    tour1[k1+1] = city21
    deleteat!(tour1, k1+2)
    tour2[k2] = city11
    tour2[k2+1] = city12
    insert!(tour2, k2+2, city13)

    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N7rr(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 2
        return Chrm
    end
    k1 = rand(1:length(tour1)-2)
    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    nt2 = length(tour2)
    Candidates = Int[] 
    if nt2 == 2  
        if city13 in Close_nodes[n_nodes+1,:] || city11 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 3  
        if city13 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        if city13 in Close_nodes[tour2[1],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, 2)
        end
    else
        if city13 in Close_nodes[n_nodes+1,:] || tour2[3] in Close_nodes[city11,:] 
            push!(Candidates, 1)
        end
        for i=2:nt2-2
            if tour2[i-1] in Close_nodes[city13,:] || tour2[i+2] in Close_nodes[city11,:]
                push!(Candidates, i)
            end
        end
        if city13 in Close_nodes[tour2[nt2-2],:] || city11 in Close_nodes[n_nodes+1,:] 
            push!(Candidates, nt2-1)
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k2 = Candidates[rand(1:length(Candidates))]
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_reverse_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
    end
    tour1[k1] = city22
    tour1[k1+1] = city21
    deleteat!(tour1, k1+2)
    tour2[k2] = city13
    tour2[k2+1] = city12
    insert!(tour2, k2+2, city11)

    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function N_cross(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Cross Exchange
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    t1 = Chrm.tours[r1].Sequence
    t2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(t1) < 6 || length(t2) < 6
        return Chrm
    end
    tau = min(8, Int(round(min(length(t1), length(t2))/2)))
    k11 = rand(1:length(t1)-tau)
    l1 = rand(1:tau)
    k12 = k11 + l1

    Candidates = Int[] 
    for i=1:length(t2)-tau
        if k11 == 1  
            if t2[i] in Close_nodes[n_nodes+1,:] 
                push!(Candidates, 1)
            end
        else
            if t2[i] in Close_nodes[t1[k11-1],:]
                push!(Candidates, i) 
            end
        end
    end

    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    k21 = Candidates[rand(1:length(Candidates))]
    l2 = rand(1:tau)
    k22 = k21 + l2
    
    new_cost1, new_cost2, straight1, straight2 = Calculate_new_cost_cross(t1, cost1, t2, cost2, k11, k12, k21, k22, T, n_nodes)
    if new_cost1 >= cost1 || new_cost2 >= cost1
        return Chrm
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

function N8(Chrm::Chromosome, TT::Matrix{Float64}, n_nodes::Int)   #divides the tours into two chuncks and attaches them
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
    if length(tour1) < 3 || length(tour2) < 3
        return Chrm
    end
    k1 = rand(2:length(tour1)-1)
    k2 = rand(2:length(tour2)-1)
    chunk11 = tour1[1:k1]
    chunk12 = tour1[k1+1:length(tour1)]
    chunk21 = tour2[1:k2]
    chunk22 = tour2[k2+1:length(tour2)]
#     chuck = 1
#     direction = true
#     distance = Inf
#     if T[chunk11[k1]+1, chunk21[1]+1] < distance 
#         distance = T[chunk11[k1]+1, chunk21[k1]+1]
#     end
#     if T[chunk11[k1]+1, chunk21[length(chunk21)]+1] < distance 
#         distance = T[chunk11[k1]+1, chunk21[length(chunk21)]+1]
#         direction = false
#     end
#     if T[chunk11[k1]+1, chunk22[1]+1] < distance 
#         distance = T[chunk11[k1]+1, chunk21[length(chunk21)]+1]
#         direction = false
#     end
    t1 = Int[]
    t2 = Int[]
    if rand() < 0.5
        if rand() < 0.5
            t1 = vcat(chunk11, chunk21)
        else
            t1 = vcat(chunk11, reverse(chunk21))
        end
        if rand() < 0.5
            t2 = vcat(chunk12, chunk22)
        else
            t2 = vcat(chunk12, reverse(chunk22))
        end
    else
        if rand() < 0.5
            t1 = vcat(chunk11, chunk22)
        else
            t1 = vcat(chunk11, reverse(chunk22))
        end
        if rand() < 0.5
            t2 = vcat(chunk12, chunk21)
        else
            t2 = vcat(chunk12, reverse(chunk21))
        end
    end
    new_cost1 = find_tour_length(t1, TT)
    if new_cost1 > cost1
        return Chrm
    end
    new_cost2 = find_tour_length(t2, TT)
    if new_cost2 > cost1
        return Chrm
    end
#     println("nice")
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.tours[r1].Sequence = t1
    Chrm.tours[r2].Sequence = t2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function Improve_chromosome(chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int, roullet::Vector{Int}, old_best::Float64)
#     Search_methods = [N1, N2, N3, N4, Ni1, Ni2, Ni3, Ni4, Ni5, Ni6, Ni7, N3r, N4sr, N4rs, N4rr, N5, N5r, N6, N6sr, N6rs, N6rr, N7, N7rs, N7sr, N7rr]
    Search_methods = [N1, Ni1, Ni2, Ni3, Ni4, Ni5]    #Ni4 not great
#     Search_methods = [N1, N2, N3, Ni1, Ni2, Ni5, Ni7, Ni3, Ni4]

    for i=1:100
        r = sample(1:length(Search_methods), weights(roullet))
#         r= rand(1:length(Search_methods))
        search = Search_methods[r]
        f1 = chrm.fitness
        chrm = search(chrm, TT, Close_nodes, demands, W, n_nodes)
#         if round(chrm.fitness, digits=4) < round(old_best, digits=4)
#             println("Improvement by local search: " , r ,"  ", round(old_best, digits=4) ," to ", round(chrm.fitness,digits=4))
#             old_best = chrm.fitness
#         end
        if chrm.fitness < f1
            roullet[r] +=1
            return chrm, r
        end
    end

    return chrm, 0
end

function Improve_Population(P::Vector{Chromosome}, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)
#     Search_methods = [N4, N2, N3, Ni6, Ni7, N3r, N4sr, N4rs, N4rr, N5, N5r, N6, N6sr, N6rs, N6rr, N7, N7rs, N7sr, N7rr]
    Search_methods = [N1, N2, N3, N4, Ni1, Ni2, Ni3, Ni4, Ni5, Ni6, Ni7, N3r, N4sr, N4rs, N4rr, N5, N5r, N6, N6sr, N6rs, N6rr, N7, N7rs, N7sr, N7rr, N_cross]
#     seq = sort(sample(1:length(P), 10, replace = false))[1:5]
#     seq = [1]
#     for j in seq
#         chrm = P[j]
#     for chrm in P[1:10]
#         for i=1:5000
#             r= rand(1:length(Search_methods))
#             search = Search_methods[r]
#             chrm = search(chrm, TT, Close_nodes, demands, W, n_nodes)
#         end

#     end
    chrm = P[1]
    P[1] = N_cross_extensive(chrm, T, n_nodes, 5)
end
        
