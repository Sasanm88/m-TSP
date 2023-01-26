using JuMP, Gurobi


function SPLIT_MIP(TT::Matrix{Float64}, demands::Vector{Int}, K::Int, W::Int, S::Vector{Int})
    
    model = Model(Gurobi.Optimizer)
    N = size(TT)[1]-2 

    @variable(model, x[1:N, 1:N] >= 0, Bin)
    @variable(model, C_max >=0)

    @objective(model, Min, C_max)
    for i=1:N
        @constraint(model, C_max >= (TT[1,S[i]+1]+TT[S[i]+1,N+2])*x[i,i])
    end

    for i=1:N-1
        for j=i+1:N
            @constraint(model, C_max >= (TT[1,S[i]+1] + sum(TT[S[k]+1, S[k+1]+1] for k=i:j-1)+ TT[S[j]+1,N+2]) * x[i,j])
        end
    end

    for i=2:N
        for j=1:i-1
            @constraint(model, x[i,j] == 0)
        end
    end
    
    @constraint(model, sum(x[1,j] for j=1:N)==1)
    @constraint(model, sum(x[i,N] for i=1:N)==1)
    
    for k=1:N-1
        @constraint(model, sum(x[i,k] for i=1:k) == sum(x[k+1,j] for j=k+1:N))
    end
    
    @constraint(model, sum(x[i,j] for i=1:N,j=1:N) <= K)
    
    for i = 1:N
        for j=i:N
            @constraint(model, sum(demands[S[k]] for k=i:j) <= W)
        end
    end

    run_time_limit = 60 #sec
    set_time_limit_sec(model, run_time_limit)
    optimize!(model)
    
    tours = Vector{Vector{Int}}()
    i = 1
    X = value.(x)
    while (i < n)
        tour = Int[]
        for j=i:n
            if X[i,j] > 0
                for t=i:j
                    push!(tour, t)
                end
                i = j+1
                break
            end
        end
        push!(tours, tour)
    end
    return objective_value(model), tours
end
