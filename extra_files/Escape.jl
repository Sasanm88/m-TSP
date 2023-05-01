
function Es1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r2].cost
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

    insert!(tour2, k2, city1)
    new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, k1)
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end


function Es2(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(1,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
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

function Es3(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,2)
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

    new_cost2 = Calculate_new_cost_add_two(tour2, cost2, city1, city2, k2, TT, n_nodes)

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

function Es4(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,2)
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

    new_cost2 = Calculate_new_cost_add_two(tour2, cost2, city2, city1, k2, TT, n_nodes)
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


function Es5(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
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


function Es6(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
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

function Es7(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
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
    city11 = tour1[k1]
    city12 = tour1[k1+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two_straight_reverse(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)

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

function Es8(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(2,2)
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

function Es9(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,3)
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

    new_cost2 = Calculate_new_cost_add_three(tour2, cost2, city1, city2, city3, k2, TT, n_nodes)

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

function Es10(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,3)
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

    new_cost2 = Calculate_new_cost_add_three(tour2, cost2, city3, city2, city1, k2, TT, n_nodes)

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

function Es11(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
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

    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)

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

function Es12(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
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

    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]

    new_cost1, new_cost2 = Calculate_new_cost_swap_three_reverse_straight(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)

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

function Es13(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
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

    city11 = tour1[k1]
    city12 = tour1[k1+1]
    city13 = tour1[k1+2]

    new_cost1, new_cost2 = Calculate_new_cost_swap_three_straight_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)

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

function Es14(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,3)
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

    city21 = tour2[k2]
    city22 = tour2[k2+1]
    city23 = tour2[k2+2]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_reverse_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22, city23, k2, TT, n_nodes)

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

function Es15(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
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
    if nt2 <= 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  
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
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)

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

function Es16(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
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
    if nt2 <= 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  
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
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_reverse_straight(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)

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

function Es17(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
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
    if nt2 <= 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  
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
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_straight_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)

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

function Es18(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(3,2)
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
    if nt2 <= 2  
        if city11 in Close_nodes[n_nodes+1,:] || city13 in Close_nodes[n_nodes+1,:]
            push!(Candidates, 1)
        end
    elseif nt2 == 4  
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
    city21 = tour2[k2]
    city22 = tour2[k2+1]
    
    new_cost1, new_cost2 = Calculate_new_cost_swap_three_with_two_reverse_reverse(tour1, cost1, city11, city12, city13, k1, tour2, cost2, city21, city22,  k2, TT, n_nodes)

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



function Es19(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Reinsert
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end

    tour1 = Chrm.tours[r1].Sequence
    cost1 = Chrm.tours[r1].cost

    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    nt = length(tour1)
    if nt <= 1
        return Chrm
    end
    Candidates = Int[] 
    if nt == 2
        Candidates = [1,2]
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour1[1] in Close_nodes[city1,:] 
            push!(Candidates, 1)
        end
        for i=2:nt-1
            if tour1[i-1] in Close_nodes[city1,:] || tour1[i+1] in Close_nodes[city1,:]
                push!(Candidates, i)
            end
        end
        if city1 in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[city1,:] 
            push!(Candidates, nt)
        end
    end
    Candidates = collect(setdiff(Set(Candidates),Set([k1])))
    if length(Candidates) == 0
        return Chrm
    end
    
    k2 = Candidates[rand(1:length(Candidates))]

    new_cost1 = Calculate_new_cost_exchange_one(tour1, cost1, city1, k1, k2, TT, n_nodes)

    deleteat!(tour1, k1)
    insert!(tour1, k2, city1)

    Chrm.tours[r1].cost = new_cost1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end


function Es20(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Exchange (permutation between two customers)
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 1
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1 = rand(1:nt)
    city1 = tour1[k1]

    Candidates = Int[] 

    for i=1:nt
        if i != k1
            if i == 1
                if city1 in Close_nodes[n_nodes+1,:] || tour1[2] in Close_nodes[city1,:]
                    push!(Candidates, 1)
                end
            elseif i == nt
                if city1 in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[city1,:]
                    push!(Candidates, nt)
                end
            else
                if tour1[i-1] in Close_nodes[city1,:] || tour1[i+1] in Close_nodes[city1,:]
                    push!(Candidates, i)
                end
            end
        end
    end
    if length(Candidates) == 0
        return Chrm
    end
    k2 = Candidates[rand(1:length(Candidates))]

    city2 = tour1[k2]

    new_cost1= Calculate_new_cost_exchange_two(tour1, cost1, city1, k1, city2, k2, TT, n_nodes)

    tour1[k1] = city2
    tour1[k2] = city1
    Chrm.tours[r1].cost = new_cost1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    return Chrm
end

function Es21(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Or-opt2 
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
    k1 = rand(1:nt-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]

    Candidates = Int[] 
    if nt == 3
        Candidates = [1,2]
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour1[1] in Close_nodes[city2,:] 
            push!(Candidates, 1)
        end
        for i=2:nt-2
            if tour1[i-1] in Close_nodes[city1,:] || tour1[i+1] in Close_nodes[city2,:]
                push!(Candidates, i)
            end
        end
        if city2 in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[city1,:] 
            push!(Candidates, nt-1)
        end
    end
    Candidates = collect(setdiff(Set(Candidates),Set([k1])))
    if length(Candidates) == 0
        return Chrm
    end

    k2 = Candidates[rand(1:length(Candidates))]

    z1 = Calculate_new_cost_or_opt2(tour1, cost1, city1, k1, city2, k2, T, n_nodes)

    deleteat!(tour1, [k1,k1+1])
    insert!(tour1, k2, city1)
    insert!(tour1, k2+1, city2)
    Chrm.tours[r1].cost = z1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end


function Es22(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Or-opt3 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 3
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1 = rand(1:nt-2)
    city1 = tour1[k1]
    city2 = tour1[k1+1]
    city3 = tour1[k1+2]
    Candidates = Int[] 
    if nt == 4
        Candidates = [1,2]
    else
        if city1 in Close_nodes[n_nodes+1,:] || tour1[1] in Close_nodes[city3,:] 
            push!(Candidates, 1)
        end
        for i=2:nt-3
            if tour1[i-1] in Close_nodes[city1,:] || tour1[i+2] in Close_nodes[city3,:]
                push!(Candidates, i)
            end
        end
        if city3 in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[city1,:] 
            push!(Candidates, nt-2)
        end
    end
    Candidates = collect(setdiff(Set(Candidates),Set([k1])))
    if length(Candidates) == 0
        return Chrm
    end

    k2 = Candidates[rand(1:length(Candidates))]
    new_cost1 = Calculate_new_cost_or_opt3(tour1, cost1, city1, city2, city3, k1, k2, T, n_nodes)

    deleteat!(tour1, [k1, k1+1, k1+2])
    insert!(tour1, k2, city1)
    insert!(tour1, k2+1, city2)
    insert!(tour1, k2+2, city3)
    Chrm.tours[r1].cost = new_cost1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Es23(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #2-opt 
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
    # indices = sample(1:length(tour1), 2, replace = false)
    i1 = rand(1:length(tour1))
    
    Candidates = Int[] 
    
    if i1 == 1
        for i=2:nt-1
            if tour1[i] in Close_nodes[n_nodes+1,:] || tour1[1] in Close_nodes[tour1[i1+1],:] 
                push!(Candidates, i)
            end
        end
    elseif i1 == nt 
        for i=2:nt-1
            if tour1[i] in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[tour1[i1-1],:] 
                push!(Candidates, i)
            end
        end
    else
        if tour1[i1] in Close_nodes[n_nodes+1,:] || tour1[1] in Close_nodes[tour1[i1+1],:] 
            push!(Candidates, 1)
        end
        if tour1[i1] in Close_nodes[n_nodes+1,:] || tour1[nt] in Close_nodes[tour1[i1-1],:] 
            push!(Candidates, nt)
        end
        for i=2:nt-1
            if i > i1 
                if tour1[i1] in Close_nodes[tour1[i+1],:] || tour1[i] in Close_nodes[tour1[i1-1],:]
                    push!(Candidates, i)
                end
            elseif i < i1
                if tour1[i1] in Close_nodes[tour1[i-1],:] || tour1[i] in Close_nodes[tour1[i1+1],:]
                    push!(Candidates, i)
                end
            end
        end
    end

    if length(Candidates) == 0
        return Chrm
    end
    Candidates = collect(Set(Candidates))
    i2 = Candidates[rand(1:length(Candidates))]
    k1, k2 = min(i1, i2), max(i1,i2)
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

function Es24(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #3-opt 
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
    k1, k2, k3 = sort!(sample(1:length(tour1), 3, replace = false))

    new_cost = Calculate_new_cost_3_opt(tour1, cost1, k1, k2, k3, T, n_nodes)

    if k2 - k1 >=3
        tour1[k1+1:k2-1] = reverse(tour1[k1+1:k2-1])
    end
    if k3 - k2 >=3
        tour1[k2+1:k3-1] = reverse(tour1[k2+1:k3-1])
    end
    Chrm.tours[r1].cost = new_cost
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Es25(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #3-permute 
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
    k1 = rand(1:length(tour1)-2)
    temp1 = copy(tour1[k1:k1+2])
    temp2 = shuffle(temp1)
    
    new_cost = Calculate_new_cost_3_permute(tour1, cost1, temp1, temp2, k1, T, n_nodes)

    tour1[k1:k1+2] = temp2
    Chrm.tours[r1].cost = new_cost
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Escape_local_optima(P::Vector{Chromosome}, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int, turn::Int64, iter::Int, max_size::Int, allowed_diff::Float64, elite::Int)
    
    local_chrm = deepcopy(P[1])
    
    methods = [Es1, Es2, Es3, Es4, Es5, Es6, Es7, Es8, Es9, Es10, Es11, Es12, Es13, Es14, Es15, Es16, Es17, Es18, Es19, Es20, Es21, Es22, Es23, Es24, Es25]
    substitues = Vector{Chromosome}()
    fs = Vector{Vector{Float64}}()

    push!(substitues, local_chrm)
    push!(fs, sort(round.([tour.cost for tour in local_chrm.tours], digits=4)))

    best_f = local_chrm.fitness

    for i=1:iter*turn
        chrm = substitues[rand(1:length(substitues))]
        chrm1 = deepcopy(chrm)
        r = 0
        r = rand(1:length(methods))
        search = methods[r]
        chrm1 = search(chrm1, TT, Close_nodes, demands, W, n_nodes)
        f = chrm1.fitness 
        if f < best_f
            if length(substitues) < max_size
                pushfirst!(substitues, chrm1)
                pushfirst!(fs, sort(round.([tour.cost for tour in chrm1.tours], digits=4)))
            else
                pushfirst!(substitues, chrm1)
                pushfirst!(fs, sort(round.([tour.cost for tour in chrm1.tours], digits=4)))
                pop!(substitues)
                pop!(fs)
            end
            best_f = f
        elseif f > best_f
            f1 = sort(round.([tour.cost for tour in chrm1.tours], digits=4))
            if (f-best_f)/best_f < allowed_diff 
                if !(f1 in fs)
                    if length(substitues) < max_size
                        inserted = false
                        for j = 1:min(elite, length(substitues))
                            if f < substitues[j].fitness
                                insert!(substitues, j, chrm1)
                                insert!(fs, j, f1)
                                inserted = true
                                break
                            end
                        end
                        if !inserted
                            push!(substitues, chrm1)
                            push!(fs, f1)
                        end
                    else
                        inserted = false
                        for j = 1:elite
                            if f < substitues[j].fitness
                                insert!(substitues, j, chrm1)
                                insert!(fs, j, f1)
                                inserted = true
                                break
                            end
                        end
                        rr = rand(elite+1:length(substitues))
                        if inserted
                            deleteat!(substitues, rr)
                            deleteat!(fs, rr)
                        else
                            substitues[rr] = chrm1
                            fs[rr] = f1
                        end
                    end
                end
            end
      
        end
    end
    sort!(substitues, by=x->x.fitness)
    for i=1:length(substitues)
        chrm = substitues[i]
        if chrm.fitness < local_chrm.fitness
            push!(P, chrm)
        else
            break
        end
    end
    sort!(P, by=x -> x.fitness)
    return substitues, fs
end

        
