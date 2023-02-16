function Creat_Random_Cromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i = 1:n_nodes])
    chromosome
end

function Generate_initial_population(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int)
    Population = Chromosome[]
    n_nodes = length(demands)
    for i=1:mu
        S = Creat_Random_Cromosome(n_nodes)
        obj, trips = SPLIT(TT, demands, K, W, S)
        push!(Population, Chromosome(S, obj, 0.0, trips))
    end
    sort!(Population, by=x -> x.fitness)
    return Population, Population[1].fitness
end

function Diversify(Population::Vector{Chromosome}, TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int)
    n_nodes = length(demands)
    n_best = Int(round(0.3 * mu))
    for i=n_best+1:length(Population)
        S = Creat_Random_Cromosome(n_nodes)
        obj, trips = SPLIT(TT, demands, K, W, S)
        Population[i] = Chromosome(S, obj, 0.0, trips)
    end
    sort!(Population, by=x -> x.fitness)
end


function Find_Closeness(TT::Matrix{Float64}, h::Float64)
    n_nodes = size(TT)[1] - 2
    num = Int(ceil(h * n_nodes))
    ClosenessT = zeros(Int, n_nodes+1, num)
    @inbounds for i = 2:n_nodes+2
        a = copy(TT[i, 2:n_nodes+1])
        b = sortperm(a)
        ClosenessT[i-1, :] = b[2:num+1]
    end
    return ClosenessT
end