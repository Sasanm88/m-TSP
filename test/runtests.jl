include(joinpath(@__DIR__, "../src/main.jl"))
    

# Solving the instances of Set 1 (Tabel 2 in the paper)
num_customers = 50   # This is N in table 2
num_salesmen = 5     # This is m in table 2
num_instances = 100  # This is the number of instances you want to solve for each combination of N and m (in our paper, we solved 100 instances)
verbose = false      # Set this to true if you want to print the detailed results for each instance
solve_instances_set1(num_customers, num_salesmen, num_instances, verbose)


# Solving the instances of Set 2 (Tabel 3 in the paper)
instances = [:eil51]  # can be set to [:eil51, :berlin52, :eil76, :rat99] for solving all instances at once
Ms = [2]              # this is the number of salemen (m in table 3) - set to [2,3,5,7] if you want to solve all the instance at once
verbose = false       # Set this to true if you want to print the detailed results for each instance
solve_instances_set2(instances, Ms, verbose)


# Solving the instances of Set 3 (Tabel 4 in the paper)
dir_name = "set3"              # This directory includes all the instances of Set 3
sample_names = ["mtsp51_3"]   # This can include as many as instances to solve in one run, e.g. ["mtsp51_3", "mtsp51_5, "mtsp100_20"] 
verbose = false                 # Set this to true if you want to print the detailed results for each instance

solve_instances_set3_4(dir_name, sample_names, verbose)

# Solving the instances of Set 4 (Tabel 5 in the paper)
dir_name = "set4"             # This directory includes all the instances of Set 4
sample_names = ["gtsp150_3"]   # This can include as many as instances to solve in one run, e.g. ["kroa200_10", "lin318_3"] 
verbose = false                # Set this to true if you want to print the detailed results for each instance

solve_instances_set3_4(dir_name, sample_names, verbose)

