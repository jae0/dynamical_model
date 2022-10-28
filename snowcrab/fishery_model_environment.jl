

if false
  # if doing manual startup .. this is done automatically on start but in case it fails:
  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
end


 
# translate model-specific functions, etc to generics
if model_variation=="logistic_discrete_basic"
  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
end


if model_variation=="logistic_discrete" 
  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
end


if model_variation=="logistic_discrete_map" 
  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
end


if model_variation=="size_structured_dde" 
  fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )  
end

if model_variation=="size_structured_dde_unnormalized" 
  fn_env = joinpath( project_directory, "size_structured_dde_unnormalized_environment.jl" )  
end

if model_variation=="size_structured_dde_ratios" 
  fn_env = joinpath( project_directory, "size_structured_dde_ratios_environment.jl" )  
end

include(  fn_env )



