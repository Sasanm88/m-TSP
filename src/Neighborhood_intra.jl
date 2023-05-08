

function Ni1!(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #Reinsert
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end

    tour1 = Chrm.tours[r1].Sequence
    cost1 = Chrm.tours[r1].cost

    
    nt = length(tour1)
    if nt <= 1
        return
    end
    
    k1 = rand(1:length(tour1))
    city1 = tour1[k1]
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
        return 
    end
    
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:length(tour1))
    new_cost1 = calculate_new_cost_exchange_one(tour1, cost1, city1, k1, k2, TT, n_nodes)
    
    if new_cost1 >= cost1
        return 
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
end


function Ni2!(Chrm::Chromosome, TT::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #Exchange (permutation between two customers)
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 1
        return 
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
        return
    end
    k2 = Candidates[rand(1:length(Candidates))]
#     k2 = rand(1:nt)
    city2 = tour1[k2]

    new_cost1= calculate_new_cost_exchange_two(tour1, cost1, city1, k1, city2, k2, TT, n_nodes)
    if new_cost1 >= cost1 
        return 
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
end

function Ni3!(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #Or-opt2 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return
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
        return
    end
    
#     k2 = rand(1:length(tour1)-1)   #Way to improve 
    k2 = Candidates[rand(1:length(Candidates))]

    z1 = calculate_new_cost_or_opt2(tour1, cost1, city1, k1, city2, k2, T, n_nodes)

    if z1 >= cost1 
        return
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

    return
end


function Ni4!(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #Or-opt3 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 3
        return
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
        return
    end

    k2 = Candidates[rand(1:length(Candidates))]
    new_cost1 = calculate_new_cost_or_opt3(tour1, cost1, city1, city2, city3, k1, k2, T, n_nodes)

    if new_cost1 >= cost1 
        return
    end
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

    return
end

function Ni5!(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #2-opt 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return
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
        return
    end
    Candidates = collect(Set(Candidates))
    i2 = Candidates[rand(1:length(Candidates))]
    k1, k2 = min(i1, i2), max(i1,i2)
    new_cost = calculate_new_cost_2_opt(tour1, cost1, k1, k2, T, n_nodes)

    if new_cost >= cost1 
        return
    end
    tour1[k1:k2] = reverse(tour1[k1:k2])
    Chrm.tours[r1].cost = new_cost
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

    return
end

function Ni6!(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #3-opt 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1, k2, k3 = sort!(sample(1:length(tour1), 3, replace = false))

    new_cost = calculate_new_cost_3_opt(tour1, cost1, k1, k2, k3, T, n_nodes)
    if new_cost >= cost1 
        return
    end
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
end

function Ni7!(Chrm::Chromosome, T::Matrix{Float64}, Close_nodes::Matrix{Int}, n_nodes::Int)   #3-permute 
    r1 = 1
    if rand() < 0.5
        r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    else
        r1 = rand(1:length(Chrm.tours))
    end
    tour1 = Chrm.tours[r1].Sequence
    if length(tour1) <= 2
        return
    end
    cost1 = Chrm.tours[r1].cost
    nt = length(tour1)
    k1 = rand(1:length(tour1)-2)
    temp1 = copy(tour1[k1:k1+2])
    temp2 = shuffle(temp1)
    
    new_cost = calculate_new_cost_3_permute(tour1, cost1, temp1, temp2, k1, T, n_nodes)
    
    if new_cost >= cost1 
        return
    end
    tour1[k1:k1+2] = temp2
    Chrm.tours[r1].cost = new_cost
    Chrm.genes = Int[]
    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    for tour in Chrm.tours
        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
    end

end