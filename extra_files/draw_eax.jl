using Plots

function Draw_tour(p1::Vector{Int}, depot::Vector{Float64}, Customers::Matrix{Float64})
    labels = Vector{Vector{String}}()

    x1 = [depot[1]]
    y1 = [depot[2]]
    labels1 = ["0"]
    for i in p1
        push!(x1, Customers[i,1])
        push!(y1, Customers[i,2])
        push!(labels1, string(i))
    end
    push!(x1, depot[1])
    push!(y1, depot[2])

    min_x = minimum(x1)
    min_y = minimum(y1)

    max_x = maximum(x1)
    max_y = maximum(y1)

    magnify = 0.2
    min_x = min_x - magnify*(max_x-min_x)
    max_x = max_x + magnify*(max_x-min_x)
    min_y = min_y - magnify*(max_y-min_y)
    max_y = max_y + magnify*(max_y-min_y)
    tit = "tsp"
    p = plot(x1, y1, marker =:circle, title = tit, label = "tour", xlim=(min_x, max_x), ylim=(min_y, max_y))

    annotate!.(x1[1:length(x1)-1]*1.02, y1[1:length(y1)-1]*1.02, text.(labels1, :left,8))
    return p
end

function Draw_nodes(p1::Vector{Int}, depot::Vector{Float64}, Customers::Matrix{Float64})
    labels = Vector{Vector{String}}()

    x1 = [depot[1]]
    y1 = [depot[2]]
    labels1 = ["0"]
    for i in 1:size(Customers)[1]
        push!(x1, Customers[i,1])
        push!(y1, Customers[i,2])
        push!(labels1, string(i))
    end
    push!(x1, depot[1])
    push!(y1, depot[2])

    min_x = minimum(x1)
    min_y = minimum(y1)

    max_x = maximum(x1)
    max_y = maximum(y1)

    magnify = 0.2
    min_x = min_x - magnify*(max_x-min_x)
    max_x = max_x + magnify*(max_x-min_x)
    min_y = min_y - magnify*(max_y-min_y)
    max_y = max_y + magnify*(max_y-min_y)
    tit = "Nodes"
    p = scatter(x1, y1, marker =:circle, title = tit, label = "tour", xlim=(min_x, max_x), ylim=(min_y, max_y))

    annotate!.(x1[1:length(x1)-1]*1.02, y1[1:length(y1)-1]*1.02, text.(labels1, :left,8))
    return p
end

# using PyPlot

function Draw_edges(E_Set::Vector{Tuple{Int64, Int64}}, Nodes::Matrix{Float64})
    x = Vector{Vector{Float64}}()
    y = Vector{Vector{Float64}}()
    min_x = Inf
    min_y = Inf
    max_x = 0.0
    max_y = 0.0
    for edge in E_Set
        x1, y1 = Nodes[edge[1]+1,:]
        x2, y2 = Nodes[edge[2]+1,:]
        push!(x, [x1, x2])
        push!(y, [y1, y2])
        if min_x > min(x1, x2)
            min_x = min(x1,x2)
        end
        if min_y > min(y1, y2)
            min_y = min(y1,y2)
        end
        if max_x < max(x1, x2)
            max_x = max(x1,x2)
        end
        if max_y < max(y1, y2)
            max_y = max(y1,y2)
        end
    end

    magnify = 0.2
    min_x = min_x - magnify*(max_x-min_x)
    max_x = max_x + magnify*(max_x-min_x)
    min_y = min_y - magnify*(max_y-min_y)
    max_y = max_y + magnify*(max_y-min_y)
    tit = "intermediate solution"
    p = plot()
    for i=1:length(x)
        p = plot!(x[i], y[i], marker =:circle, title = tit, label="", xlim=(min_x, max_x), ylim=(min_y, max_y))
    end

    labels = [string(i) for i=0:size(Nodes)[1]-1]
    annotate!.(Nodes[:,1]*1.02, Nodes[:,2]*1.02, text.(labels, :left,8))
    return p
end

function Draw_cycles(cycles::Vector{Vector{Int}}, Nodes::Matrix{Float64})

    x = Vector{Vector{Float64}}()
    y = Vector{Vector{Float64}}()
    labels = Vector{Vector{String}}()

    min_x = Inf
    max_x = 0.0
    min_y = Inf
    max_y = 0.0
    for cycle in cycles
        x1 = Float64[]
        y1 = Float64[]
        labels1 = String[]
        for i in cycle
            push!(x1, Nodes[i+1,1])
            push!(y1, Nodes[i+1,2])
            push!(labels1, string(i))
        end
        push!(x1, x1[1])
        push!(y1, y1[1])
        push!(x, x1)
        push!(y, y1)
        push!(labels, labels1)
        if minimum(x1) < min_x
            min_x = minimum(x1)
        end
        if minimum(y1) < min_y
            min_y = minimum(y1)
        end
        if maximum(x1) > max_x
            max_x = maximum(x1)
        end
        if maximum(y1) > max_y
            max_y = maximum(y1)
        end
    end
    magnify = 0.2
    min_x = min_x - magnify*(max_x-min_x)
    max_x = max_x + magnify*(max_x-min_x)
    min_y = min_y - magnify*(max_y-min_y)
    max_y = max_y + magnify*(max_y-min_y)
    tit = "Cycels"
    p = plot(x[1], y[1], marker =:circle, title = tit, label = "c1", xlim=(min_x, max_x), ylim=(min_y, max_y))
    for i=1:length(x)
        annotate!.(x[i][1:length(x[i])-1], y[i][1:length(y[i])-1], text.(labels[i], :left,8))
    end
    for i=2:length(x)
        p = plot!(x[i], y[i], marker =:circle, label = "c"*string(i))
    end 
    return p
end