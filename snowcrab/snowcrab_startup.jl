
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
    outputs_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "output", "fishery_model" ) 
  end

  mkpath(outputs_directory)
#   cd( outputs_directory )   # this is necessary as julia stores packages (versions) specific to this project here 

  print( "outputs_directory: ", outputs_directory )
  mkpath(model_outdir)
  print( "outputs_directory: ", outputs_directory )


# ---------------
# make a copy of the input data in case ... 

  if  occursin( r"size_structured", model_variation ) 
    fndat_source = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", 
      "1999_present_fb", "fishery_model_results", "turing1", "biodyn_number_size_struct.RData" )
  elseif  occursin( r"logistic_discrete", model_variation ) 
    fndat_source = joinpath( homedir(), "bio.data", "bio.snowcrab", "modelled", 
      "1999_present_fb", "fishery_model_results", "turing1", "biodyn_biomass.RData" )
  end

  fndat = joinpath( model_outdir, basename(fndat_source) )
  cp( fndat_source, fndat; force=true )

  
 
  if  occursin( r"logistic_discrete", model_variation ) 
        pkgs = [
            "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random",
            "ForwardDiff", "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",  
            "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
            "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
        ]
  end

    # add Turing@v0.21.10  # to add a particular version
    #  "DynamicHMC", 

  if occursin( r"size_structured", model_variation ) 
        pkgs = [
            "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "QuadGK",
            "MCMCChains", "DynamicPPL", "AdvancedHMC", "DistributionsAD", "Bijectors",  
            "AbstractPPL", "Memoization", # "Enzyme", "Diffractor",
            "ForwardDiff", "DataFrames", "CSV", "JLD2", "PlotThemes", "Colors", "ColorSchemes", "RData", 
            "MKL", "Plots", "StatsPlots",  "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
            "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
        ]
   end


# ---------------
# LOAD environment (libs and functions)
  if  occursin( r"size_structured", model_variation ) 
    fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )
  elseif  occursin( r"logistic_discrete", model_variation ) 
    fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  end
  


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
