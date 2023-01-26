using Distances
using Random


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
