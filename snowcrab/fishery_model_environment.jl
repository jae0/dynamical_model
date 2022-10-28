

if false
  # if doing manual startup .. this is done automatically on start but in case it fails:
  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
end


 
# translate model-specific functions, etc to generics
 
if  occursin( r"size_structured", model_variation ) 

  fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )
  
elseif  occursin( r"logistic_discrete", model_variation ) 
  
  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  

end

include(  fn_env )



