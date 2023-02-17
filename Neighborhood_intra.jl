function Calculate_new_cost_exchange_one(tour::Vector{Int}, cost::Float64, city::Int, position1::Int, 
    position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    if position1 == 1
        if position2 == nt 
            new_cost = cost - T[1, city+1] - T[city+1, tour[2]+1] + T[1, tour[2]+1] + T[city+1, n_nodes+2] 
            + T[tour[nt]+1, city+1] - T[tour[nt]+1, n_nodes+2]
        else
            new_cost = cost - T[1, city+1] - T[city+1, tour[2]+1] + T[1, tour[2]+1] + T[tour[position2]+1, city+1] 
            + T[city+1, tour[position2+1]+1] - T[tour[position2]+1, tour[position2+1]+1] 
        end
    elseif position1 == nt
        if position2 == 1 
            new_cost = cost + T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1] - T[city+1, n_nodes+2] 
            - T[tour[nt-1]+1, city+1] + T[tour[nt-1]+1, n_nodes+2]
        else
            new_cost = cost + T[tour[position2-1]+1, city+1] + T[city+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1] 
            - T[city+1, n_nodes+2] - T[tour[nt-1]+1, city+1] + T[tour[nt-1]+1, n_nodes+2]
        end
    else
        if position2 == nt 
            new_cost = cost - T[tour[position1-1]+1, city+1] - T[city+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
             + T[city+1, n_nodes+2] + T[tour[nt]+1, city+1] - T[tour[nt]+1, n_nodes+2]
        else
            new_cost = cost - T[tour[position1-1]+1, city+1] - T[city+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
            + T[tour[position2]+1, city+1] + T[city+1, tour[position2+1]+1] - T[tour[position2]+1, tour[position2+1]+1] 
        end
    end
    return new_cost
end

function Ni1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])

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
    if new_cost1 >= cost1
        return Chrm
    end

    deleteat!(tour1, k1)
    insert!(tour1, k2, city1)

    Chrm.tours[r1].cost = new_cost1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
    # print("Ni1  ")
    return Chrm
end