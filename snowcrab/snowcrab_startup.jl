

  if ! @isdefined project_directory 
    # defaulting to location of this file "bio.snowcrab/inst/julia" 
    # my call is: JULIA_NUM_THREADS=4 julia -i ~/projects/dynamical_model/snowcrab/startup.jl
    project_directory = @__DIR__() 
  end
 
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  import Pkg  # or using Pkg
  Pkg.activate(project_directory)  # so now you activate the package
  Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use
  print( "project_directory: ", project_directory )


  import Pkg  # or using Pkg
 
  Pkg.activate(project_directory)  # so now you activate the package

  Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use
 
  if ! @isdefined outputs_directory 
    # tailor to your specific installation
    outputs_directory = joinpath( bio_data_directory, "bio.snowcrab",  "fishery_model" ) 
  end

  mkpath(outputs_directory)
#   cd( outputs_directory )   # this is necessary as julia stores packages (versions) specific to this project here 
  print( "outputs_directory: ", outputs_directory )

  model_outdir = joinpath( outputs_directory, string(year_assessment), model_variation )
  mkpath(model_outdir)

  print( "model_outdir: ", model_outdir )
 
 
  if  occursin( r"logistic_discrete", model_variation ) 
        pkgs = [
            "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "Setfield", "Memoization",
            "ForwardDiff", "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",  
            "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
            "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
        ]
  end 

    # add Turing@v0.21.10  # to add a particular version
    #  "DynamicHMC", 

  if occursin( r"size_structured", model_variation ) 
        pkgs = [
            "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "QuadGK", "Setfield", "Memoization",
            "MCMCChains", "DynamicPPL", "AdvancedHMC", "DistributionsAD", "Bijectors",  
            "AbstractPPL", "Memoization", # "Enzyme", "Diffractor",
            "ForwardDiff", "DataFrames", "CSV", "JLD2", "PlotThemes", "Colors", "ColorSchemes", "RData", 
            "MKL", "Plots", "StatsPlots",  "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
            "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
        ]
  end

  # load libs and check settings
  # pkgs are defined in snowcrab_startup.jl
  for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages

# ---------------
# LOAD environment (libs and functions)
  if  occursin( r"size_structured", model_variation ) 
    fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )
  elseif  occursin( r"logistic_discrete", model_variation ) 
    fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  end
  

  include( fn_env )  # loads libs and setup workspace / data (fn_env is defined in the snowcrab_startup.jl)
    
#=

project_directory = @__DIR__() #  same folder as the file
push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg

Pkg.activate(project_directory)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

cd( project_directory )

print( project_directory )

=#
