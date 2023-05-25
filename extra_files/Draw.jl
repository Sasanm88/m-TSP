using Plots


function Draw_Solution(chrm::Chromosome, depot::Vector{Float64}, Customers::Matrix{Float64}, obj::Float64)
    x = Vector{Vector{Float64}}()
    y = Vector{Vector{Float64}}()
    labels = Vector{Vector{String}}()
    m = length(chrm.tours)
    min_x = Inf
    max_x = 0.0
    min_y = Inf
    max_y = 0.0
    for j=1:m 
        seq = chrm.tours[j].sequence
        x1 = [depot[1]]
        y1 = [depot[2]]
        labels1 = ["0"]
        for i in seq
            push!(x1, Customers[i,1])
            push!(y1, Customers[i,2])
            push!(labels1, string(i))
        end
        push!(x1, depot[1])
        push!(y1, depot[2])
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
    tit = "minmax objective = " * string(round(obj, digits=1))
    p = plot(x[1], y[1], marker =:circle, title = tit, label = "tour1", xlim=(min_x, max_x), ylim=(min_y, max_y))
    for i=1:m
        annotate!.(x[i][1:length(x[i])-1].+20, y[i][1:length(y[i])-1].+20, text.(labels[i], :left,6))
    end
    for i=2:m
        p = plot!(x[i], y[i], marker =:circle, label = "tour"*string(i))
    end 
    for (j,tour) in enumerate(chrm.tours)
        println("Tour ", j, ":" , tour.cost)
#         for (k,i) in enumerate(tour.Sequence)
#             print(i, " ")
#             if k%25==0
#                 println()
#             end
#         end
#         println()
#         print("cost=", tour.cost)
#         println()
    end
    return p
end