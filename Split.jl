mutable struct Label
    Ri::Vector{Int}
    Vir::Vector{Float64}
    Pir::Vector{Int}
    Cir::Vector{Float64}
end


function SPLIT(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, S::Vector{Int}) #In m-TSP, demands is a vector of ones and W is infinity
    n = length(demands)
    Labels = Label[]
    for i=1:n
        R = Int[]
        if i==n
            R = [i for i=1:K]
        else
            R = [i for i=1:min(i,K-1)]
        end
        V = Float64[]
        P = Int[]
        C = Float64[]
        for r in R
            push!(V, Inf)
            push!(P, n+1)
            push!(C, Inf)
        end
        push!(Labels, Label(R, V, P, C))
    end


    for i=1:n
        R = [0]
        if i > 1
            R = Labels[i-1].Ri
        end
        for r in R
            Current_V = 0
            if i > 1
                Current_V = Labels[i].Vir[r]
            end
            if Current_V < Inf    #This is V^i_r
                load = 0
                t = 0
                j = i
                while (j<=n) && (load<=W)
                    load += demands[S[j]]
                    if i==j
                        t = TT[1,S[j]+1] + TT[S[j]+1, n+2]
                    else
                        t = t - TT[S[j-1]+1,1] + TT[S[j-1]+1, S[j]+1] + TT[S[j]+1, n+2]
                    end
                    if load<=W && r+1 in Labels[j].Ri
                        old_t = 0.0
                        if i > 1
                            old_t = Labels[i-1].Vir[r]
                        end
                        new_t = max(old_t, t)
                        if new_t < Labels[j].Vir[r+1]
                            Labels[j].Vir[r+1] = new_t
                            Labels[j].Pir[r+1] = i-1
                            Labels[j].Cir[r+1] = t
                        end
                    end
                    j += 1
                end
            end
        end
    end

    trips = Vector{Tour}()

    rs = argmin(Labels[n].Vir)
    for i = 1:rs
        push!(trips, Tour(Int[], 0.0))
    end
    t = rs
    j = n
    while t>0
        i = Labels[j].Pir[t]
        trips[t].cost = Labels[j].Cir[t] 
        for k=i+1:j
            push!(trips[t].Sequence, S[k])
        end
        t -= 1
        j = i
    end
    obj = minimum(Labels[n].Vir)
    return obj, trips
end


