function Calculate_new_cost_exchange_one(tour::Vector{Int}, cost::Float64, city::Int, position1::Int, 
    position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    
    t1 = copy(tour)
    deleteat!(t1, position1)
    insert!(t1, position2, city)
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    
    
#     if position1 == 1
#         if position2 == nt 
#             new_cost = cost - T[1, city+1] - T[city+1, tour[2]+1] + T[1, tour[2]+1] + T[city+1, n_nodes+2] 
#             + T[tour[nt]+1, city+1] - T[tour[nt]+1, n_nodes+2]
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("A")
#             end
#         else
#             new_cost = cost - T[1, city+1] - T[city+1, tour[2]+1] + T[1, tour[2]+1] + T[tour[position2]+1, city+1] 
#             + T[city+1, tour[position2+1]+1] - T[tour[position2]+1, tour[position2+1]+1] 
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("B")
#             end
#         end
#     elseif position1 == nt
#         if position2 == 1 
#             new_cost = cost + T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1] - T[city+1, n_nodes+2] 
#             - T[tour[nt-1]+1, city+1] + T[tour[nt-1]+1, n_nodes+2]
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("C")
#             end
#         else
#             new_cost = cost + T[tour[position2-1]+1, city+1] + T[city+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1] 
#             - T[city+1, n_nodes+2] - T[tour[nt-1]+1, city+1] + T[tour[nt-1]+1, n_nodes+2]
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("D")
#             end
#         end
#     else
#         if position2 == nt 
#             new_cost = cost - T[tour[position1-1]+1, city+1] - T[city+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
#              + T[city+1, n_nodes+2] + T[tour[nt]+1, city+1] - T[tour[nt]+1, n_nodes+2]
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("E")
#             end
#         else
#             new_cost = cost - T[tour[position1-1]+1, city+1] - T[city+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
#             + T[tour[position2-1]+1, city+1] + T[city+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1] 
#             if round(new_cost, digits=4) != round(z1, digits=4)
#                 println("F")
#             end
#         end
#     end
    return z1
end

function Calculate_new_cost_exchange_two(tour::Vector{Int}, cost::Float64, city1::Int, position1::Int, city2::Int,
    position2::Int, T::Matrix{Float64}, n_nodes::Int)
    
    t1 = copy(tour)
    t1[position1] = city2
    t1[position2] = city1
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    
#     nt = length(tour)
#     if position1 == 1
#         if position2 == nt 
#             new_cost = cost - T[1, city1+1] - T[city1+1, tour[2]+1] + T[1, city2+1] + T[city2+1, tour[2]+1] 
#             + T[city1+1, n_nodes+2] + T[tour[nt-1]+1, city1+1] - T[city2+1, n_nodes+2] - T[tour[nt-1]+1, city2+1]
#         else
#             new_cost = cost - T[1, city1+1] - T[city1+1, tour[2]+1] + T[1, city2+1] + T[city2+1, tour[2]+1]
#             + T[city1+1, tour[position2+1]+1] + T[tour[position2-1]+1, city1+1] - T[city2+1, tour[position2+1]+1] - T[tour[position2-1]+1, city2+1]
#         end
#     elseif position1 == nt
#         if position2 == 1 
#             new_cost = cost - T[1, city2+1] - T[city2+1, tour[2]+1] + T[1, city1+1] + T[city1+1, tour[2]+1] 
#             + T[city2+1, n_nodes+2] + T[tour[nt-1]+1, city2+1] - T[city1+1, n_nodes+2] - T[tour[nt-1]+1, city1+1]
#         else
#             new_cost = cost - T[tour[nt-1]+1, city1+1] - T[city1+1, n_nodes+2] + T[tour[nt-1]+1, city2+1] + T[city2+1, n_nodes+2]
#             + T[city1+1, tour[position2+1]+1] + T[tour[position2-1]+1, city1+1] - T[city2+1, tour[position2+1]+1] - T[tour[position2-1]+1, city2+1]
#         end
#     else
#         if position2 == 1 
#             new_cost = cost - T[1, city2+1] - T[city2+1, tour[2]+1] + T[1, city1+1] + T[city1+1, tour[2]+1] 
#             + T[city2+1, tour[position1+1]+1] + T[tour[position1-1]+1, city2+1] - T[city1+1, tour[position1+1]+1] - T[tour[position1-1]+1, city1+1]
#         elseif position2 == nt 
#             new_cost = cost - T[tour[nt-1]+1, city2+1] - T[city2+1, n_nodes+2] + T[tour[nt-1]+1, city1+1] + T[city1+1, n_nodes+2]
#             + T[city2+1, tour[position1+1]+1] + T[tour[position1-1]+1, city2+1] - T[city1+1, tour[position1+1]+1] - T[tour[position1-1]+1, city1+1]
#         else
#             new_cost = cost + T[city1+1, tour[position2+1]+1] + T[tour[position2-1]+1, city1+1] - T[city2+1, tour[position2+1]+1] - T[tour[position2-1]+1, city2+1]
#             + T[city2+1, tour[position1+1]+1] + T[tour[position1-1]+1, city2+1] - T[city1+1, tour[position1+1]+1] - T[tour[position1-1]+1, city1+1]
#         end
#     end
    return z1
end

function Ni1(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Reinsert
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
#     k2 = rand(1:length(tour1))
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


function Ni2(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Exchange (permutation between two customers)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
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
#     k2 = rand(1:nt)
    city2 = tour1[k2]

    new_cost1= Calculate_new_cost_exchange_two(tour1, cost1, city1, k1, city2, k2, TT, n_nodes)
    if new_cost1 >= cost1 
        return Chrm
    end
    tour1[k1] = city2
    tour1[k2] = city1
    Chrm.tours[r1].cost = new_cost1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end
#     print("Ni2 ")
    return Chrm
end

function Ni3(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Or-opt2 
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1 = rand(1:nt-1)
    city1 = tour1[k1]
    city2 = tour1[k1+1]

    k2 = rand(1:length(tour1)-1)

    t1 = copy(tour1)
    deleteat!(t1, [k1,k1+1])
    insert!(t1, k2, city1)
    insert!(t1, k2+1, city2)
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    if z1 >= cost1 
        return Chrm
    end
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


function Ni4(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #Or-opt3 
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
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
    k2 = rand(1:length(tour1)-2)

    t1 = copy(tour1)
    deleteat!(t1, [k1, k1+1, k1+2])
    insert!(t1, k2, city1)
    insert!(t1, k2+1, city2)
    insert!(t1, k2+2, city3)
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    if z1 >= cost1 
        return Chrm
    end
    deleteat!(tour1, [k1, k1+1, k1+2])
    insert!(tour1, k2, city1)
    insert!(tour1, k2+1, city2)
    insert!(tour1, k2+2, city3)
    Chrm.tours[r1].cost = z1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Ni5(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #2-opt 
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    indices = sample(1:length(tour1), 2, replace = false)
    k1, k2 = minimum(indices), maximum(indices)

    t1 = copy(tour1)
    t1[k1:k2] = reverse(t1[k1:k2])
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    if z1 >= cost1 
        return Chrm
    end
    tour1[k1:k2] = reverse(tour1[k1:k2])
    Chrm.tours[r1].cost = z1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Ni6(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #3-opt 
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1, k2, k3 = sort!(sample(1:length(tour1), 3, replace = false))

    t1 = copy(tour1)
    if k2 - k1 >=3
        t1[k1+1:k2-1] = reverse(t1[k1+1:k2-1])
    end
    if k3 - k2 >=3
        t1[k2+1:k3-1] = reverse(t1[k2+1:k3-1])
    end
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    if z1 >= cost1 
        return Chrm
    end
    if k2 - k1 >=3
        tour1[k1+1:k2-1] = reverse(tour1[k1+1:k2-1])
    end
    if k3 - k2 >=3
        tour1[k2+1:k3-1] = reverse(tour1[k2+1:k3-1])
    end
    Chrm.tours[r1].cost = z1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end

function Ni7(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, demands::Vector{Int}, W::Int, n_nodes::Int)   #3-permute 
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return Chrm
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1= rand(1:length(tour1)-2)

    t1 = copy(tour1)
    t1[k1:k1+2] = shuffle!(t1[k1:k1+2])
    temp = copy(t1[k1:k1+2])    
    pushfirst!(t1, 0)
    push!(t1, n_nodes+1)
    z1 = 0.0
    for i=1:length(t1)-1
        z1 += T[t1[i]+1, t1[i+1]+1]
    end
    if z1 >= cost1 
        return Chrm
    end
    tour1[k1:k1+2] = temp
    Chrm.tours[r1].cost = z1
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return Chrm
end