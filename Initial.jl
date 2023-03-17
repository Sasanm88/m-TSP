using LKH

function find_tsp_tour1(Ct::Matrix{Float64})
    scale_factor = 1000
    dist_mtx = round.(Int, Ct .* scale_factor)

    # tsp_tour, tsp_tour_len = Concorde.solve_tsp(dist_mtx)
    tsp_tour, tour_length = LKH.solve_tsp(dist_mtx)

    @assert tsp_tour[1] == 1

    return tsp_tour[2:length(tsp_tour)].-1, tour_length/scale_factor
end

function Change_initial(c::Vector{Int64}, n_nodes::Int64)
    cc = copy(c)
    r = rand()
    if r < 0.5
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
        if idx1 > idx2
            temp = idx1
            idx1 = idx2
            idx2 = temp
        end
        cc[idx1:idx2] = reverse(cc[idx1:idx2])
    else
        idx = sample(1:n_nodes, rand(2:Int(floor(n_nodes/2))), replace = false) 
        cc[idx] = shuffle(cc[idx])
    end
    return cc
end

function Creat_Random_Cromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i = 1:n_nodes])
    chromosome
end

function Generate_initial_population(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int, tsp_tour::Vector{Int})
    Population = Chromosome[]
    n_nodes = length(demands)
    obj, trips = SPLIT(TT, demands, K, W, tsp_tour)
    push!(Population, Chromosome(tsp_tour, obj, 0.0, trips))
    S = Int[]
    for i=1:mu-1
        if rand() < 0.5
            S = Change_initial(tsp_tour, n_nodes)
        else
            S = Creat_Random_Cromosome(n_nodes)
        end
        obj, trips = SPLIT(TT, demands, K, W, S)
        push!(Population, Chromosome(S, obj, 0.0, trips))
    end
    sort!(Population, by=x -> x.fitness)
    return Population, Population[1].fitness
end

function Diversify(Population::Vector{Chromosome}, TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, mu::Int, tsp_tour::Vector{Int})
    n_nodes = length(demands)
    n_best = Int(round(0.15 * mu)) 
    for i=n_best+1:length(Population)
        S = Int[]
        if rand() < 0.5
            chrm = deepcopy(Population[1])
            ch = N8(chrm, TT, n_nodes)
            Population[i] = ch
        else
            S = Creat_Random_Cromosome(n_nodes)
            obj, trips = SPLIT(TT, demands, K, W, S)
            Population[i] = Chromosome(S, obj, 0.0, trips)
        end
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