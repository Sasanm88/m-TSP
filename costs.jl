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
            new_cost1 += T[tour1[position1-1]+1, city21+1] + T[city21+1, city22+1] + T[city22+1, tour1[position1+2]+1] - T[tour1[position1-1]+1, city11+1] - T[city11+1, city12+1] - T[city12+1, tour1[position1+2]+1]
        end
    end
    
#     t2 = copy(tour2)
#     t2[position2] = city11
#     t2[position2+1] = city12
#     pushfirst!(t2, 0)
#     push!(t2, n_nodes+1)

#     z2 = 0.0
#     for i=1:length(t2)-1
#         z2 += T[t2[i]+1, t2[i+1]+1]
#     end
    
    if nt2 == 2
        new_cost2 = T[1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2] 
    else
        if position2 == 1
            new_cost2 += T[1, city11+1] + T[city11+1, city12+1]  + T[city12+1, tour2[3]+1] - T[1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[3]+1]
        elseif position2 == nt2-1
            new_cost2 += T[tour2[nt2-2]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, n_nodes+2] - T[tour2[nt2-2]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, n_nodes+2]
        else
            new_cost2 += T[tour2[position2-1]+1, city11+1] + T[city11+1, city12+1] + T[city12+1, tour2[position2+2]+1] - T[tour2[position2-1]+1, city21+1] - T[city21+1, city22+1] - T[city22+1, tour2[position2+2]+1]
        end
    end
    return new_cost1, new_cost2
end


function Calculate_new_cost_exchange_one(tour1::Vector{Int}, cost::Float64, city::Int, position1::Int, 
    position2::Int, T::Matrix{Float64}, n_nodes::Int)
    tour = copy(tour1)
    nt = length(tour)
    if position1 == position2
        return cost
    end
    if position1 == 1
        cost = cost - T[1, city+1] - T[city+1, tour[2]+1] + T[1, tour[2]+1]
    elseif position1 == nt
        cost = cost - T[city+1, n_nodes+2] - T[tour[nt-1]+1, city+1] + T[tour[nt-1]+1, n_nodes+2]
    else
        cost = cost - T[tour[position1-1]+1, city+1] - T[city+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
    end

    deleteat!(tour, position1)
    if position2 == 1
        cost = cost + T[1, city+1] + T[city+1, tour[1]+1] - T[1, tour[1]+1]
    elseif position2 == nt
        cost = cost + T[tour[nt-1]+1, city+1] + T[city+1, n_nodes+2] - T[tour[nt-1]+1, n_nodes+2]
    else
        cost = cost + T[tour[position2-1]+1, city+1] + T[city+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1]
    end
    return cost
end

function Calculate_new_cost_exchange_two(tour::Vector{Int}, cost::Float64, city1::Int, position1::Int, city2::Int,
    position2::Int, T::Matrix{Float64}, n_nodes::Int)
    nt = length(tour)
    
    if position1 == position2
        return cost
    end
    if nt==2
        cost = T[1, tour[2]+1] + T[tour[2]+1, tour[1]+1] + T[tour[1]+1, n_nodes+2]
        return cost
    end
    if position2 == position1 +1
        if position1 == 1
            cost = cost - T[1, city1+1] - T[city2+1, tour[3]+1] + T[1,city2+1] + T[city1+1, tour[3]+1]
        elseif position2 == nt
            cost = cost - T[city2+1, n_nodes+2] - T[tour[nt-2]+1, city1+1] + T[city1+1, n_nodes+2] + T[tour[nt-2]+1, city2+1]
        else
            cost = cost - T[tour[position1-1]+1, city1+1] - T[city2+1, tour[position2+1]+1] + T[tour[position1-1]+1, city2+1] + T[city1+1, tour[position2+1]+1] 
        end
        return cost
    end
    if position1 == position2 +1
        if position2 == 1
            cost = cost - T[1, city2+1] - T[city1+1, tour[3]+1] + T[1,city1+1] + T[city2+1, tour[3]+1]
        elseif position1 == nt
            cost = cost - T[city1+1, n_nodes+2] - T[tour[nt-2]+1, city2+1] + T[city2+1, n_nodes+2] + T[tour[nt-2]+1, city1+1]
        else
            cost = cost - T[tour[position2-1]+1, city2+1] - T[city1+1, tour[position1+1]+1] + T[tour[position2-1]+1, city1+1] + T[city2+1, tour[position1+1]+1] 
        end
        return cost
    end
    if position1 == 1
        cost = cost - T[1, city1+1] - T[city1+1, tour[2]+1] + T[1, tour[2]+1]
        cost = cost + T[1, city2+1] + T[city2+1, tour[2]+1] - T[1, tour[2]+1]
    elseif position1 == nt
        cost = cost - T[city1+1, n_nodes+2] - T[tour[nt-1]+1, city1+1] + T[tour[nt-1]+1, n_nodes+2]
        cost = cost + T[tour[nt-1]+1, city2+1] + T[city2+1, n_nodes+2] - T[tour[nt-1]+1, n_nodes+2]
    else
        cost = cost + T[tour[position1-1]+1, city2+1] + T[city2+1, tour[position1+1]+1] - T[tour[position1-1]+1, tour[position1+1]+1]
        cost = cost - T[tour[position1-1]+1, city1+1] - T[city1+1, tour[position1+1]+1] + T[tour[position1-1]+1, tour[position1+1]+1]
    end

    if position2 == 1
        cost = cost - T[1, city2+1] - T[city2+1, tour[2]+1] + T[1, tour[2]+1]
        cost = cost + T[1, city1+1] + T[city1+1, tour[2]+1] - T[1, tour[2]+1]
    elseif position2 == nt
        cost = cost - T[city2+1, n_nodes+2] - T[tour[nt-1]+1, city2+1] + T[tour[nt-1]+1, n_nodes+2]
        cost = cost + T[tour[nt-1]+1, city1+1] + T[city1+1, n_nodes+2] - T[tour[nt-1]+1, n_nodes+2]
    else
        cost = cost + T[tour[position2-1]+1, city1+1] + T[city1+1, tour[position2+1]+1] - T[tour[position2-1]+1, tour[position2+1]+1]
        cost = cost - T[tour[position2-1]+1, city2+1] - T[city2+1, tour[position2+1]+1] + T[tour[position2-1]+1, tour[position2+1]+1]
    end
    return cost
end

function Calculate_new_cost_or_opt2(tour1::Vector{Int}, cost::Float64, city1::Int, position1::Int, 
    city2::Int, position2::Int, T::Matrix{Float64}, n_nodes::Int)
    
    if position1 == position2
        return cost
    end
    tour = copy(tour1)
    nt = length(tour)
    if position1 == 1
        cost = cost - T[1, city1+1] - T[city2+1, tour[3]+1] + T[1, tour[3]+1]
    elseif position1 == nt-1
        cost = cost - T[city2+1, n_nodes+2] - T[tour[nt-2]+1, city1+1] + T[tour[nt-2]+1, n_nodes+2]
    else
        cost = cost - T[tour[position1-1]+1, city1+1] - T[city2+1, tour[position1+2]+1] + T[tour[position1-1]+1, tour[position1+2]+1]
    end

    deleteat!(tour, [position1, position1+1])
    
    if position2 == 1
        cost = cost + T[1, city1+1] + T[city2+1, tour[1]+1] - T[1, tour[1]+1]
    elseif position2 == nt-1
        cost = cost + T[tour[nt-2]+1, city1+1] + T[city2+1, n_nodes+2] - T[tour[nt-2]+1, n_nodes+2]
    else
        cost = cost + T[tour[position2-1]+1, city1+1] + T[city2+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1]
    end
    return cost
end


function Calculate_new_cost_or_opt3(tour1::Vector{Int}, cost::Float64, city1::Int, city2::Int, city3::Int, position1::Int, 
     position2::Int, T::Matrix{Float64}, n_nodes::Int)
    
    if position1 == position2
        return cost
    end
    tour = copy(tour1)
    nt = length(tour)
    if position1 == 1
        cost = cost - T[1, city1+1] - T[city3+1, tour[4]+1] + T[1, tour[4]+1]
    elseif position1 == nt-2
        cost = cost - T[city3+1, n_nodes+2] - T[tour[nt-3]+1, city1+1] + T[tour[nt-3]+1, n_nodes+2]
    else
        cost = cost - T[tour[position1-1]+1, city1+1] - T[city3+1, tour[position1+3]+1] + T[tour[position1-1]+1, tour[position1+3]+1]
    end

    deleteat!(tour, [position1, position1+1, position1+2])
    
    if position2 == 1
        cost = cost + T[1, city1+1] + T[city3+1, tour[1]+1] - T[1, tour[1]+1]
    elseif position2 == nt-2
        cost = cost + T[tour[nt-3]+1, city1+1] + T[city3+1, n_nodes+2] - T[tour[nt-3]+1, n_nodes+2]
    else
        cost = cost + T[tour[position2-1]+1, city1+1] + T[city3+1, tour[position2]+1] - T[tour[position2-1]+1, tour[position2]+1]
    end
    return cost
end

function Calculate_new_cost_2_opt(tour::Vector{Int}, cost::Float64, position1::Int,
        position2::Int, T::Matrix{Float64}, n_nodes::Int)

    nt = length(tour)
    if position1 == 1
        cost = cost - T[1, tour[1]+1] + T[1,tour[position2]+1]
    else
        cost = cost - T[tour[position1-1]+1, tour[position1]+1] + T[tour[position1-1]+1, tour[position2]+1]
    end
    
    if position2 == nt
        cost = cost - T[tour[nt]+1, n_nodes+2] + T[tour[position1]+1, n_nodes+2]
    else
        cost = cost - T[tour[position2]+1, tour[position2+1]+1] + T[tour[position1]+1, tour[position2+1]+1]
    end
    return cost
end

function Calculate_new_cost_3_opt(tour::Vector{Int}, cost::Float64, k1::Int, k2::Int, k3::Int,
         T::Matrix{Float64}, n_nodes::Int)

    nt = length(tour)
    if k3 - k2 > 2
        position1 = k2+1
        position2 = k3-1
        cost = cost - T[tour[position1-1]+1, tour[position1]+1] + T[tour[position1-1]+1, tour[position2]+1] - T[tour[position2]+1, tour[position2+1]+1] + T[tour[position1]+1, tour[position2+1]+1]
    end
    if k2 - k1 > 2
        position1 = k1+1
        position2 = k2-1
        cost = cost - T[tour[position1-1]+1, tour[position1]+1] + T[tour[position1-1]+1, tour[position2]+1] - T[tour[position2]+1, tour[position2+1]+1] + T[tour[position1]+1, tour[position2+1]+1]
    end
    return cost
end

function Calculate_new_cost_3_permute(tour::Vector{Int}, cost::Float64, S1::Vector{Int}, S2::Vector{Int}
        ,k1::Int, T::Matrix{Float64}, n_nodes::Int)
    if S1 == S2 
        return cost
    end
    nt = length(tour)
    if nt == 3
        cost = T[1, S2[1]+1] + T[S2[1]+1, S2[2]+1] + T[S2[2]+1, S2[3]+1] + T[S2[3]+1, n_nodes+2]
        return cost
    end
    cost = cost - T[S1[1]+1, S1[2]+1] - T[S1[2]+1, S1[3]+1] + T[S2[1]+1, S2[2]+1] + T[S2[2]+1, S2[3]+1]
    
    if k1 == 1
        cost = cost - T[1, S1[1]+1] + T[1, S2[1]+1] - T[S1[3]+1, tour[4]+1] + T[S2[3]+1, tour[4]+1] 
    elseif k1 == nt-2
        cost = cost - T[S1[3]+1, n_nodes+2] + T[S2[3]+1, n_nodes+2] - T[tour[nt-3]+1, S1[1]+1] + T[tour[nt-3]+1, S2[1]+1]
    else
        cost = cost - T[tour[k1-1]+1, S1[1]+1] + T[tour[k1-1]+1, S2[1]+1] - T[S1[3]+1, tour[k1+3]+1] + T[S2[3]+1, tour[k1+3]+1] 
    end
    return cost
end