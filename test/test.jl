include(joinpath(@__DIR__, "../src/main.jl")) 

instances = [:eil51]
Ms = [2]
test(instances, Ms)

dir_name = "set1"
sample_names = ["mtsp150_3", "mtsp150_5", "mtsp150_10", "kroa200_3", "kroa200_5","kroa200_10", "lin318_3", "lin318_5", "lin318_10"]

Solve_instances(dir_name, sample_names)