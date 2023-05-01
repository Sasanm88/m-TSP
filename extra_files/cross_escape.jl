function N_cross_extensive(Chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int, tau_::Int)   #Cross Exchange
    
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    t1 = Chrm.tours[r1].Sequence
    cost1 = Chrm.tours[r1].cost
    if length(t1) < 4
        return Chrm
    end
    routes = [i for i=1:length(Chrm.tours)]
    for r2 in setdiff(routes, r1)
        t2 = Chrm.tours[r2].Sequence
        cost2 = Chrm.tours[r2].cost
        tau_p = min(tau_, min(length(t1), length(t2))-2)
        if length(t2) >= 4
            for tau1=1:tau_p
                for k11 = 1:length(t1)-tau1
                    k12 = k11 + tau1
                    for tau2=1:tau_p
                        for k21 = 1:length(t2)-tau2
                            k22 = k21+tau2
                            
                            new_cost1, new_cost2, straight1, straight2 = Calculate_new_cost_cross(t1, cost1, t2, cost2, k11, k12, k21, k22, T, n_nodes)
#                             if new_cost1 < cost1 && new_cost2 < cost1
#                                 println(new_cost1, "  ", new_cost2)
#                             end
                            if new_cost1 < cost1 && new_cost2 < cost1

                                if straight2
                                    alpha1 = copy(t1[k11:k12])
                                else
                                    alpha1 = reverse(copy(t1[k11:k12]))
                                end
                                if straight1
                                    alpha2 = copy(t2[k21:k22])
                                else
                                    alpha2 = reverse(copy(t2[k21:k22]))
                                end
                                deleteat!(t1, [i for i=k11:k12])
                                for i=1:k22-k21+1
                                    insert!(t1, i+k11-1, alpha2[i])
                                end

                                deleteat!(t2, [i for i=k21:k22])
                                for i=1:k12-k11+1
                                    insert!(t2, i+k21-1, alpha1[i])
                                end

                                Chrm.tours[r1].cost = new_cost1
                                Chrm.tours[r2].cost = new_cost2
                                Chrm.genes = Int[]
                                Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
                                for tour in Chrm.tours
                                    Chrm.genes = vcat(Chrm.genes, tour.Sequence)
                                end
                                return Chrm
                            end
                        end
                    end
                end
            end
        end
    end
    return Chrm
end

function N_cross_super_extensive(Chrm::Chromosome, T::Matrix{Float64}, n_nodes::Int, tau_::Int)   #Cross Exchange
    Closenodes = Find_Closeness(T, 0.3)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    t1 = Chrm.tours[r1].Sequence
    cost1 = Chrm.tours[r1].cost
    if length(t1) < 4
        return Chrm
    end
    routes = [i for i=1:length(Chrm.tours)]
    for r2 in setdiff(routes, r1)
        t2 = Chrm.tours[r2].Sequence
        cost2 = Chrm.tours[r2].cost
        tau_p = min(tau_, min(length(t1), length(t2))-2)
        if length(t2) >= 4
            for tau1=1:tau_p
                for k11 = 1:length(t1)-tau1
                    k12 = k11 + tau1
                    for tau2=1:tau_p
                        for k21 = 1:length(t2)-tau2
                            k22 = k21+tau2
#                             positions1 = intersect(vcat([pos for pos=1:k11],[pos for pos=k21+2:length(t1)-tau2]), union(Closenodes[t2[k21]],Closenodes[t2[k22]]))
                            positions1 = vcat([pos for pos=1:k11],[pos for pos=k21+2:length(t1)-tau2])
#                             positions2 = intersect(vcat([pos for pos=1:k21],[pos for pos=k22+2:length(t2)-tau1]), union(Closenodes[t1[k11]],Closenodes[t1[k12]]))
                            positions2 = vcat([pos for pos=1:k21],[pos for pos=k22+2:length(t2)-tau1])
                            for k13 in positions1
                                for k23 in positions2
                                    new_cost1, new_cost2, straight1, straight2 = Calculate_new_cost_cross_upgraded(t1, cost1, t2, cost2, k11, k12, k21, k22, k13, k23, T, n_nodes)
#                             if new_cost1 < cost1 && new_cost2 < cost1
#                                 println(new_cost1, "  ", new_cost2)
#                             end
                                    if new_cost1 < cost1 && new_cost2 < cost1
                                        print("yes")
                                        if straight2
                                            alpha1 = copy(t1[k11:k12])
                                        else
                                            alpha1 = reverse(copy(t1[k11:k12]))
                                        end
                                        if straight1
                                            alpha2 = copy(t2[k21:k22])
                                        else
                                            alpha2 = reverse(copy(t2[k21:k22]))
                                        end
                                        deleteat!(t1, [i for i=k11:k12])
                                        if k13 > k12
                                            k13 = k13 - (k12-k11+1)
                                        end
                                        for i=1:k22-k21+1
                                            insert!(t1, i+k13-1, alpha2[i])
                                        end

                                        deleteat!(t2, [i for i=k21:k22])
                                        if k23 > k22
                                            k23 = k23 - (k22-k21+1)
                                        end
                                        for i=1:k12-k11+1
                                            insert!(t2, i+k23-1, alpha1[i])
                                        end

                                        Chrm.tours[r1].cost = new_cost1
                                        Chrm.tours[r2].cost = new_cost2
                                        Chrm.genes = Int[]
                                        Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
                                        for tour in Chrm.tours
                                            Chrm.genes = vcat(Chrm.genes, tour.Sequence)
                                        end
                                        return Chrm
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return Chrm
end

function N_shift1_extensive(Chrm::Chromosome, TT::Matrix{Float64}, n_nodes::Int)   #Shift(0,1)
    r1 = argmax([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
    routes = [i for i=1:length(Chrm.tours)]
    tour1 = Chrm.tours[r1].Sequence
    cost1 = Chrm.tours[r1].cost
    for r2 in setdiff(routes, r1)
    
        tour2 = Chrm.tours[r2].Sequence
        cost2 = Chrm.tours[r2].cost
        for k1=1:length(tour1)
            city1 = tour1[k1]
            for k2=1:length(tour2)
    
                new_cost2 = Calculate_new_cost_add_one(tour2, cost2, city1, k2, TT, n_nodes)
                if new_cost2 < cost1
                    insert!(tour2, k2, city1)
                    new_cost1 = Calculate_new_cost_remove_one(tour1, cost1, k1, TT, n_nodes)
                    deleteat!(tour1, k1)
                    Chrm.tours[r1].cost = new_cost1
                    Chrm.tours[r2].cost = new_cost2
                    Chrm.genes = Int[]
                    Chrm.fitness = maximum([Chrm.tours[i].cost for i=1:length(Chrm.tours)])
                    for tour in Chrm.tours
                        Chrm.genes = vcat(Chrm.genes, tour.Sequence)
                    end
                    return Chrm
                end
            end
        end
    end
    return Chrm
end