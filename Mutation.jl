function Displacement_mutation(child::Vector{Int}, n_nodes::Int)  
    idx1 = rand(2:n_nodes-1)
    idx2 = rand(2:n_nodes-1)
    while idx1 == idx2
        idx1 = rand(2:n_nodes-1)
        idx2 = rand(2:n_nodes-1)
    end

    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end

    subtour = child[idx1:idx2]
    deleteat!(child, idx1:idx2)
    idx3 = rand(1:length(child))

    for i = 1:idx2-idx1+1
        insert!(child, idx3 + i, subtour[i])
    end
end

function Exchange_mutation(child::Vector{Int}, n_nodes::Int)
    idx1, idx2 = sample(1:n_nodes, Weights(ones(n_nodes)), 2, replace=false)
    child[idx1], child[idx2] = child[idx2], child[idx1]
end


function Insertion_mutation(child::Vector{Int}, n_nodes::Int)
    idx1 = rand(1:n_nodes)
    city = child[idx1]
    deleteat!(child, idx1)
    idx2 = rand(1:length(child))
    insert!(child, idx2, city)
end

function Simple_Inversion_mutation(child::Vector{Int}, n_nodes::Int)
    idx1 = rand(1:n_nodes)
    idx2 = rand(1:n_nodes)
    while idx1 == idx2
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
    end

    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    child[idx1:idx2] = reverse(child[idx1:idx2])
end


function Inversion_mutation(child::Vector{Int}, n_nodes::Int)
    idx1 = rand(2:n_nodes-1)
    idx2 = rand(2:n_nodes-1)
    while idx1 == idx2
        idx1 = rand(2:n_nodes-1)
        idx2 = rand(2:n_nodes-1)
    end

    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end

    subtour = child[idx1:idx2]
    deleteat!(child, idx1:idx2)
    idx3 = rand(1:length(child))
    reverse!(subtour)
    for i = 1:idx2-idx1+1
        insert!(child, idx3 + i, subtour[i])
    end
end


function Scramble_mutation(child::Vector{Int}, n_nodes::Int)
    idx1 = rand(1:n_nodes)
    idx2 = rand(1:n_nodes)
    while idx1 == idx2
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
    end

    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    child[idx1:idx2] = shuffle!(child[idx1:idx2])
end

function tour_mutation(child::Vector{Int}, n_nodes::Int)
    idx = sample(1:n_nodes, Int(round(0.2 * n_nodes)), replace=false)
    child[idx] = child[shuffle(idx)]
end


function Mutate(child::Vector{Int}, Mutation_Chance::Float64)
    n_nodes = length(child)
    Mutation_functions = [Displacement_mutation, Exchange_mutation, Insertion_mutation, Simple_Inversion_mutation, Inversion_mutation, Scramble_mutation, tour_mutation]
    if rand() < Mutation_Chance
        Mutation_functions[rand(1:length(Mutation_functions))](child, n_nodes)
    end
end
