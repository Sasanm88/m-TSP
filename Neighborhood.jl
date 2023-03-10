
function N1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Shift(0,1)
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

    
    if new_cost2 >= cost1
        return Chrm
    end
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


function N2(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Swap(1,1)
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
#     k2 = rand(1:length(tour2)-1)
    city21 = tour2[k2]
    city22 = tour2[k2+1]

    new_cost1, new_cost2 = Calculate_new_cost_swap_two(tour1, cost1, city11, city12, k1, tour2, cost2, city21, city22, k2, TT, n_nodes)
#     t2 = copy(tour2)
#     t2[k2] = city11
#     t2[k2+1] = city12
#     pushfirst!(t2, 0)
#     push!(t2, n_nodes+1)

#     z2 = 0.0
#     for i=1:length(t2)-1
#         z2 += T[t2[i]+1, t2[i+1]+1]
#     end
#     if round(new_cost2, digits = 6) != round(z2, digits = 6)
#         print("BBB   ")
#     end
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

function Improve_chromosome(chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int, roullet::Vector{Int})
#         Search_methods = [N1, N2, N3, N4, Ni1, Ni2, Ni3, Ni4, Ni5, Ni6, Ni7]
        Search_methods = [N1, Ni1, Ni2, Ni3, Ni4, Ni5]    #Ni4 not great
        for i=1:100
    #     shuffle!(Search_methods)
    #         for search in Search_methods
            r = sample(1:length(Search_methods), weights(roullet))
            search = Search_methods[r]
            f1 = chrm.fitness
            chrm = search(chrm, TT, Close_nodes, demands, W, n_nodes)
            if chrm.fitness < f1
    #             println("AA")
                roullet[r] +=1
            end
        end
        return chrm
    end
        
