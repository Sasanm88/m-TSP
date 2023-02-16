function Calculate_new_cost_add_one(tour::Vector{Int}, cost::Float64, city::Int, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 0 
        return T[1, city+1] + T[city+1, n_modes+2]
    end
    if position == 1
        cost += T[1,city+1]+T[city+1, tour[1]+1]-T[1,tour[1]+1]
    elseif position == nt+1
        cost += T[city+1, n_nodes+2]+T[tour[nt]+1, city+1]-T[tour[nt]+1, n_nodes+2]
    else
        cost += T[tour[position-1]+1, city+1] + T[city+1, tour[position]+1]-T[tour[position-1]+1,tour[position]+1]
    end
    return cost
end

function Calculate_new_cost_add_two(tour::Vector{Int}, cost::Float64, city1::Int, city2::Int, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 0 
        return T[1, city1+1] + T[city1+1, city2+1] + T[city2+1, n_modes+2]
    end
    if position == 1
        cost += T[1, city1+1] + T[city1+1, city2+1] + T[city2+1, tour[1]+1] - T[1,tour[1]+1]
    elseif position == nt+1
        cost += T[city2+1, n_nodes+2] + T[city1+1, city2+1] +T[tour[nt]+1, city1+1] - T[tour[nt]+1, n_nodes+2]
    else
        cost += T[tour[position-1]+1, city1+1] + T[city1+1, city2+1] + T[city2+1, tour[position]+1]-T[tour[position-1]+1,tour[position]+1]
    end
    return cost
end

function Calculate_new_cost_remove_one(tour::Vector{Int}, cost::Float64, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 1
        return 0.0
    end
    if position == 1
        cost += T[1, tour[2]+1] - T[tour[1]+1, tour[2]+1] - T[1,tour[1]+1]
    elseif position == nt
        cost += T[tour[nt-1]+1, n_nodes+2] - T[tour[nt-1]+1, tour[nt]+1] - T[tour[nt]+1, n_nodes+2]
    else
        cost += T[tour[position-1]+1, tour[position+1]+1] -T[tour[position-1]+1, tour[position]+1] - T[tour[position]+1, tour[position+1]+1]
    end
    return cost
end

function Calculate_new_cost_remove_two(tour::Vector{Int}, cost::Float64, position::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if nt == 2
        return 0.0
    end
    if position == 1
        cost += T[1, tour[3]+1] - T[tour[1]+1, tour[2]+1] - T[tour[2]+1, tour[3]+1] - T[1,tour[1]+1]
    elseif position == nt-1
        cost += T[tour[nt-2]+1, n_nodes+2] - T[tour[nt-2]+1, tour[nt-1]+1] - T[tour[nt-1]+1, tour[nt]+1] - T[tour[nt]+1, n_nodes+2]
    else
        cost += T[tour[position-1]+1, tour[position+2]+1] -T[tour[position-1]+1, tour[position]+1] - T[tour[position]+1, tour[position+1]+1] -T[tour[position+1]+1, tour[position+2]+1]
    end
    return cost
end

function Calculate_new_cost_swap_one(tour1::Vector{Int}, cost1::Float64, city1::Int, position1::Int, tour2::Vector{Int}, cost2::Float64, city2::Int, position2::Int,T::Matrix{Float64}, n_nodes::Int)
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost1 = cost1
    new_cost2 = cost2
    if nt1 == 1
        new_cost1 = T[1, city2+1] + T[city2+1, n_nodes+2] 
    else
        if position1 == 1
            new_cost1 += T[1, city2+1] + T[city2+1, tour1[2]+1] - T[1, city1+1] - T[city1+1, tour1[2]+1]
        elseif position1 == nt1
            new_cost1 += T[tour1[nt1-1]+1, city2+1] + T[city2+1, n_nodes+2] - T[tour1[nt1-1]+1, city1+1] - T[city1+1, n_nodes+2]
        else
            new_cost1 += T[tour1[position1-1]+1, city2+1] + T[city2+1, tour1[position1+1]+1] - T[tour1[position1-1]+1, city1+1] - T[city1+1, tour1[position1+1]+1]
        end
    end
    if nt2 == 1
        new_cost2 = T[1, city1+1] + T[city1+1, n_nodes+2] 
    else
        if position2 == 1
            new_cost2 += T[1, city1+1] + T[city1+1, tour2[2]+1] - T[1, city2+1] - T[city2+1, tour2[2]+1]
        elseif position2 == nt2
            new_cost2 += T[tour2[nt2-1]+1, city1+1] + T[city1+1, n_nodes+2] - T[tour2[nt2-1]+1, city2+1] - T[city2+1, n_nodes+2]
        else
            new_cost2 += T[tour2[position2-1]+1, city1+1] + T[city1+1, tour2[position2+1]+1] - T[tour2[position2-1]+1, city2+1] - T[city2+1, tour2[position2+1]+1]
        end
    end
    return new_cost1, new_cost2
end

function Calculate_new_cost_swap_two(tour1::Vector{Int}, cost1::Float64, city11::Int, city12::Int, position1::Int, tour2::Vector{Int}, cost2::Float64, city21::Int, city22::Int, position2::Int,T::Matrix{Float64}, n_nodes::Int)
    nt1 = length(tour1)
    nt2 = length(tour2)
    new_cost1 = cost1
    new_cost2 = cost2
    if nt1 == 2
        new_cost1 = T[1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2] 
    else
        if position1 == 1
            new_cost1 += T[1, city21+1] + T[city21+1, city22+1]  + T[city22+1, tour1[3]+1] - T[1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[3]+1]
        elseif position1 == nt1-1
            new_cost1 += T[tour1[nt1-2]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, n_nodes+2] - T[tour1[nt1-2]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, n_nodes+2]
        else
            new_cost1 += T[tour1[position1-1]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[position1+1]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[position1+1]+1]
        end
    end
    if nt2 == 2
        new_cost2 = T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2] 
    else
        if position2 == 1
            new_cost2 += T[1, city11+1] + T[city11+1, city12+1]  + T[city12+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
        elseif position2 == nt2-1
            new_cost2 += T[tour2[nt2-2]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
        else
            new_cost2 += T[tour2[position2-1]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, tour2[position2+1]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+1]+1]
        end
    end
    return new_cost1, new_cost2
end

function N1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r1].cost
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
    candidates = intersect(tour2, Close_nodes[city1,:])
    if length(tour2) == 0
        candidates = [1]
    end
    if length(candidates) == 0
        return Chrm
    end
    city = candidates[rand(1:length(candidates))]
    k2 = findfirst(x->x==city, tour2)
    if rand()<0.5
        k2 +=1
    end
    new_cost2 = Calculate_new_cost_add_one(tour2, cost2, city1, k2, TT, n_nodes)
    if new_cost2 >= cost1
        return Chrm
    end
    insert!(tour2, city1, k2)
    new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
    deleteat!(tour1, k1)
    Chrm.tours[r1].cost = new_cost1
    Chrm.tours[r2].cost = new_cost2
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    print("N1")
    return Chrm
end


function N2(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(1,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r1].cost
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
#     print("N2")
    return Chrm
end

function N3(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,2)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    r2 = setdiff(routes, r1)[rand(1:length(Chrm.tours)-1)]
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) < 2
        return Chrm
    end
    tour2 = Chrm.tours[r2].Sequence
    cost1 = Chrm.tours[r1].cost
    cost2 = Chrm.tours[r1].cost
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
    if new_cost2 >= cost1
        return Chrm
    end
    print("A")
    insert!(tour2, city1, k2)
    insert!(tour2, city2, k2+1)
    
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
    cost2 = Chrm.tours[r1].cost
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
        if city11 in Close_nodes[nt2-2,:] || city12 in Close_nodes[n_nodes+1,:] 
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
        if city11 in Close_nodes[nt2-2,:] || city12 in Close_nodes[n_nodes+1,:] 
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
#     print("N4  ")
    return Chrm
end

function Improve_chromosome(chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)
    Search_methods = [N1, N2, N3, N4]
    shuffle!(Search_methods)
    for search in Search_methods
        # @code_warntype N5(chrm, TT, DD, ClosenessT, ClosenessD, n_nodes)
        f1 = chrm.fitness
        chrm = search(chrm, TT, Close_nodes, demands, W, n_nodes)
#         f2 = chrm.fitness
#         if f2<f1
#             print("A")
#         end
    end
    return chrm
end
        
