function segments_intersect_orient(A::AbstractVector{Float64}, B::AbstractVector{Float64}, C::AbstractVector{Float64}, D::AbstractVector{Float64})
    # Unpack the coordinates of the first segment into separate variables
    x1, y1 = A
    x2, y2 = B

    # Unpack the coordinates of the second segment into separate variables
    x3, y3 = C
    x4, y4 = D

    # Calculate the orientations of the line segments
    o1 = orientation((x1, y1), (x2, y2), (x3, y3))
    o2 = orientation((x1, y1), (x2, y2), (x4, y4))
    o3 = orientation((x3, y3), (x4, y4), (x1, y1))
    o4 = orientation((x3, y3), (x4, y4), (x2, y2))

    # Check if the line segments intersect
    if o1 != o2 && o3 != o4
        return true
    end

    # Check if the line segments are collinear and overlapping
    #     if o1 == 0 && on_segment((x1, y1), (x3, y3), (x2, y2))
    #         return true
    #     end

    #     if o2 == 0 && on_segment((x1, y1), (x4, y4), (x2, y2))
    #         return true
    #     end

    #     if o3 == 0 && on_segment((x3, y3), (x1, y1), (x4, y4))
    #         return true
    #     end

    #     if o4 == 0 && on_segment((x3, y3), (x2, y2), (x4, y4))
    #         return true
    #     end

    # If none of the above conditions are met, the line segments do not intersect
    return false
end

function orientation(p, q, r)
    # Calculate the orientation of the triplet (p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])

    if val == 0
        # Points are collinear
        return 0
    elseif val > 0
        # Clockwise orientation
        return 1
    else
        # Counterclockwise orientation
        return 2
    end
end

function on_segment(p, q, r)
    # Check if point q lies on segment pr
    if q[0] <= max(p[0], r[0]) && q[0] >= min(p[0], r[0]) &&
       q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1])
        return true
    end

    return false
end


function find_intersections(t1::Vector{Int}, t2::Vector{Int}, Customers::Matrix{Float64}, depot::Vector{Float64})
    #     intersection_points = Vector{Vector{Int}}()
    for i = 0:length(t1)
        if i == 0 || i == length(t1)
            loopset = [k for k  in 1:length(t2)-1]
        else
            loopset = [k for k = 0:length(t2)]
        end

        if i == 0
            node1 = @view depot[:, 1]  # essentially same as `node1 = depot`, but for type stability
            node2 = @view Customers[t1[i+1], :]
        elseif i == length(t1)
            node1 = @view Customers[t1[i], :]
            node2 = @view depot[:, 1]
        else
            node1 = @view Customers[t1[i], :]
            node2 = @view Customers[t1[i+1], :]
        end    

        for j in loopset
            intersected = false

            if j == 0
                node3 = @view depot[:, 1]
                node4 = @view Customers[t2[j+1], :]
            elseif j == length(t2)
                node3 = @view Customers[t2[j], :]
                node4 = @view depot[:, 1]
            else
                node3 = @view Customers[t2[j], :]
                node4 = @view Customers[t2[j+1], :]
            end

            if segments_intersect_orient(node1, node2, node3, node4)
                return i, j, true
            end
        end
    end
    return 0, 0, false
end

function solve_one_intersections(t1::Vector{Int}, t2::Vector{Int}, T::Matrix{Float64}, k1::Int, k2::Int)
    n1 = length(t1)
    n2 = length(t2)
    if k1 == 0
        a = Int[]
    else
        a = t1[1:k1]
    end
    if k2 == n2
        b = Int[]
    else
        b = t2[k2+1:n2]
    end
    if k2 == 0
        c = Int[]
    else
        c = t2[1:k2]
    end
    if k1 == n1
        d = Int[]
    else
        d = t1[k1+1:n1]
    end
    tour11 = vcat(a, b)
    tour12 = vcat(c, d)
    tour21 = vcat(a, reverse(c))
    tour22 = vcat(reverse(b), d)
    cost11 = find_tour_length(tour11, T)
    cost12 = find_tour_length(tour12, T)
    cost21 = find_tour_length(tour21, T)
    cost22 = find_tour_length(tour22, T)
    if max(cost11, cost12) < max(cost21, cost22)
        return tour11, tour12, cost11, cost12
    else
        return tour21, tour22, cost21, cost22
    end
end

function solve_all_intersections!(chrm::Chromosome, Customers::Matrix{Float64}, depot::Vector{Float64}, T::Matrix{Float64})
    m = length(chrm.tours)
    n_nodes = length(chrm.genes)
    intersected_ = true
    while intersected_
        intersected_ = false
        for i in 1:m-1
            for j = i+1:m
                tour1 = chrm.tours[i].sequence
                tour2 = chrm.tours[j].sequence
                if length(tour1) > 0 && length(tour2) > 0
                    intersected = true
                    while intersected
                        k1, k2, intersected = find_intersections(tour1, tour2, Customers, depot)
                        if intersected
                            t1, t2, c1, c2 = solve_one_intersections(tour1, tour2, T, k1, k2)
                            tour1 = copy(t1)
                            tour2 = copy(t2)
                            chrm.tours[i].cost = c1
                            chrm.tours[j].cost = c2
                            chrm.tours[i].sequence = tour1
                            chrm.tours[j].sequence = tour2
                            intersected_ = true
                        end
                    end
                end
            end
        end
        for tour in chrm.tours
            two_opt_on_route(tour, T, n_nodes)
        end
    end

    if rand() < 0.1
        improve_after_removing_intersections(chrm.tours, T, n_nodes, m, Customers, depot)
    end
    chrm.genes = Int[]
    chrm.fitness = 0.0
    for tour in chrm.tours
        chrm.genes = vcat(chrm.genes, tour.sequence)
        if chrm.fitness < tour.cost
            chrm.fitness = tour.cost
        end
    end
end

function two_opt_on_route(tour::Tour, T::Matrix{Float64}, n_nodes::Int)   #2-opt 
    tour1 = tour.sequence
    cost1 = tour.cost
    nt = length(tour1)
    improved = true
    while improved
        improved = false
        for i1  in 1:nt-1
            if improved
                break
            end
            for i2 = i1+1:nt
                new_cost = calculate_new_cost_2_opt(tour1, cost1, i1, i2, T, n_nodes)
                if round(new_cost, digits=2) < round(cost1, digits=2)
                    cost1 = new_cost
                    improved = true
                    tour1[i1:i2] = reverse(tour1[i1:i2])
                    tour.cost = new_cost
                    break
                end
            end
        end
    end
end

function improve_after_removing_intersections(tours::Vector{Tour}, T::Matrix{Float64}, n_nodes::Int, m::Int, Customers::Matrix{Float64}, depot::Vector{Float64})
    sort!(tours, by=x -> x.cost, rev=true)
    means = [mean(Customers[t1, :], dims=1)[1, :] for t1 in [tours[i].sequence for i in 1:m]]

    distances = zeros(m, m)
    for i in 1:m-1
        for j = i+1:m
            distances[i, j] = euclidean(means[i], means[j])
            distances[j, i] = distances[i, j]
        end
    end

    all_tours = [i for i in 1:m]

    improved = true

    while improved
        current_tour = 1
        improved = false
        next_tours = sortperm(distances[current_tour, :])[2:min(m, 3)]
        #         println(tours[current_tour].cost)
        for next_tour in next_tours
            #             println(tours[next_tour].cost)
            tour1 = tours[current_tour].sequence
            cost1 = tours[current_tour].cost
            tour2 = tours[next_tour].sequence
            cost2 = tours[next_tour].cost
            for i in 1:length(tour1)
                new_cost1 = calculate_new_cost_remove_one(tour1, cost1, i, T, n_nodes)
                new_cost2 = Inf
                best_position = 0
                for j in 1:length(tour2)+1
                    temp = calculate_new_cost_add_one(tour2, cost2, tour1[i], j, T, n_nodes)
                    if temp < new_cost2
                        new_cost2 = temp
                        best_position = j
                    end
                end
                #                 println(i,"  ", best_position)

                if new_cost2 < cost1
                    t2 = copy(tour2)
                    insert!(t2, best_position, tour1[i])
                    if i == 1
                        node1 = depot
                        node2 = Customers[tour1[i+1], :]
                    elseif i == length(tour1)
                        node1 = Customers[tour1[i-1], :]
                        node2 = depot
                    else
                        node1 = Customers[tour1[i-1], :]
                        node2 = Customers[tour1[i+1], :]
                    end
                    node3 = Customers[tour1[i], :]
                    intersected = false
                    for j in 1:length(t2)
                        if j == 1
                            node4 = depot
                            node5 = Customers[t2[j], :]
                        elseif j == length(t2) + 1
                            node4 = Customers[t2[length(t2)], :]
                            node5 = depot
                        else
                            node4 = Customers[t2[j-1], :]
                            node5 = Customers[t2[j], :]
                        end

                        intersect1 = segments_intersect_orient(node1, node2, node3, node4)
                        intersect2 = segments_intersect_orient(node1, node2, node3, node5)
                        if intersect1 || intersect2
                            intersected = true
                            break
                        end
                    end
                    if !intersected
                        insert!(tour2, best_position, tour1[i])
                        deleteat!(tour1, i)
                        tours[current_tour].cost = new_cost1
                        tours[next_tour].cost = new_cost2
                        sort!(tours, by=x -> x.cost, rev=true)
                        improved = true
                        break
                    end
                end
            end
        end
    end
end