
# good writeup
# https://www.multibugs.org/documentation/latest/spatial/SpatialDistributions.html#Appendix




# in REPL or VSCODE, this needs to be loaded first as "startup.jl" is skipped
# choose one:
# project_directory = joinpath(  @__DIR__(), "Documents", "GitHub", "dynamical_model", "car" ) #  same folder as the current file
# project_directory = "/home/jae/projects/dynamical_model/car"

push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
# include( "startup.jl" )


 
using Turing 
using PDMats
using LinearAlgebra 
using StatsModels
using DataFrames
using SparseArrays
using Graphs


# ---------------
# load support functions
include( joinpath( project_directory, "car_functions.jl"  ))  


# load libs and options and prepare data for diffeq/turing model and set default parameters
include( joinpath( project_directory, "scottish_lip_cancer_environment.jl"  ))  # bootstrap different project environments depending on above choices


 


# various test implementations:

m = turing_car(D, W, X, log_offset, y)        # 204 sec

m = turing_car_prec(D, W, X, log_offset, y)   # 65 sec

# Morris' "simple_iar" testing difference formulation .. 
# using same run specs: results are similar with much better ess than stan   
# m = turing_icar_direct_test( node1, node2, std(y) )  #

m = turing_icar_direct_bym(X, log_offset, y, node1, node2) # 13 sec
   
# bym2 requires a "scaling factor" 
m = turing_icar_direct_bym2(X, log_offset, y, node1, node2, scaling_factor_bym2(W)) # W is adjacency matrix , 18 sec; 30min for full run


if false
  # bym2 grouped model (multiple groups or disconnected groups): this is not finished 
  groups_unique = unique(sort(groups))
  gi = Vector{Vector{Int64}}()
  for g in groups_unique
      o =  findall(x -> x==g, groups) 
      push!(gi, o)
  end

  scaling_factor = scaling_factor_bym2_groups(node1, node2, groups)
  m = turing_icar_direct_bym2_groups(X, log_offset, y, node1, node2, scaling_factor, groups)

end





# check timings and accuracy
n_samples, n_adapts, n_chains = 100, 100, 1

# use NUTS: see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
target_acceptance, max_depth, init_ϵ = 0.65, 7, 0.05
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)

# sample
o = sample(m, turing_sampler, n_samples) 




# larger runs  
n_samples, n_adapts, n_chains = 9_000, 1_000, 4
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.1   # Morris uses 0.97 for target_acceptance, stan default is 0.95; such high acceptance rate does not work well -- divergent chains
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)
o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
# o = sample(m, turing_sampler, n_samples) 

# if on windows and threads are still not working, use single processor mode:
# o = mapreduce(c -> sample(m, turing_sampler, n_samples), chainscat, 1:n_chains)

 

showall( summarize(o) )

 

 