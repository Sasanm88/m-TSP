using Distances
using Random
using TSPLIB

function Creat_Sample(n::Int, length_of_area::Int, tspeed::Int)
    #(Number of demand nodes, length of square area (km), truck speed (km/h)
    depot = (rand(0:1000*length_of_area), rand(0:1000*length_of_area))
    Nodes = Vector{Tuple{Int, Int}}()

    for i=1:n
        push!(Nodes, (rand(0:1000*length_of_area), rand(0:1000*length_of_area)))
    end
    T = zeros(n,n)
    Tp = zeros(n)

    for i=1:n

        Tp[i] = cityblock(depot, Nodes[i])*60/(1000*tspeed)

        for j=1:n
            T[i,j] = cityblock(Nodes[i], Nodes[j])*60/(1000*tspeed)
        end
    end

    TT = zeros(n+2,n+2)
  
    TT[2:n+1,2:n+1] = T
    TT[2:n+1,1] = Tp
    TT[1,2:n+1] = Tp
    TT[2:n+1,n+2] = Tp
    TT[n+2,2:n+1] = Tp

    return TT
end

function Calculate_distance_matrices_TSPLIB(tspeed::Int, allNodes::Matrix{Float64})
    num_of_nodes = size(allNodes)[1] - 1
    T = Matrix{Float64}(undef, num_of_nodes + 2, num_of_nodes + 2)
    depot = allNodes[1, :]
    Nodes = allNodes[2:num_of_nodes+1, :]

    T[1, 1] = 0.0
    T[num_of_nodes+2, num_of_nodes+2] = 0.0
    T[1, num_of_nodes+2] = 0.0
    T[num_of_nodes+2, 1] = 0.0
    @inbounds for i in 1:num_of_nodes
        T[1, i+1] = euclidean(depot, Nodes[i, :]) / tspeed
        T[i+1, 1] = T[1, i+1]
        T[num_of_nodes+2, i+1] = T[1, i+1]
        T[i+1, num_of_nodes+2] = T[1, i+1]
        @inbounds for j in 1:num_of_nodes
            T[i+1, j+1] = euclidean(Nodes[i, :], Nodes[j, :]) / tspeed
            T[j+1, i+1] = T[i+1, j+1]
        end
    end

    return T
end

function Read_TSPLIB_instance(sample_name::Symbol, tspeed::Int)
    tsp = readTSPLIB(sample_name)
    dEligible = Int[]
    T = Calculate_distance_matrices_TSPLIB(tspeed, tsp.nodes)
    return T
end