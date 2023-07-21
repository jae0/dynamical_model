
  print( "\n\nWARNING: if this is the initial run, it will take a while to precompile/download all libraries. \n\n" )

  if ! @isdefined project_directory 
    # defaulting to location of this file "bio.snowcrab/inst/julia" 
    # my call is: JULIA_NUM_THREADS=4 julia -i ~/projects/dynamical_model/snowcrab/startup.jl
    project_directory = @__DIR__() 
  end
 
  print( "project_directory: ", project_directory, "\n\n" )


  import Pkg  # or using Pkg
  Pkg.activate(project_directory)  # so now you activate the package
  Base.active_project()  
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  cd( project_directory )

  
   
  if ! @isdefined outputs_directory 
    # tailor to your specific installation
    outputs_directory = joinpath( bio_data_directory, "bio.snowcrab",  "fishery_model" ) 
  end

  mkpath(outputs_directory)
  print( "outputs_directory: ", outputs_directory, "\n\n" )

  model_outdir = joinpath( outputs_directory, string(year_assessment), model_variation )
  mkpath(model_outdir)

  print( "model_outdir: ", model_outdir, "\n\n" )
 
 
  if  occursin( r"logistic_discrete", model_variation ) 
        pkgs = [
            "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random", "Setfield", "Memoization",
            "ForwardDiff", "DataFrames", "JLD2", "CSV", "PlotThemes", "Colors", "ColorSchemes", "RData",  
            "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
            "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
        ]
  end  

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
 
  # for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages

  include( "startup.jl" )


  function discretize_decimal( x, delta=0.01 ) 
    num_digits = Int(ceil( log10(1.0 / delta)) )   # time floating point rounding
    out = round.( round.( x ./ delta; digits=0 ) .* delta; digits=num_digits)
    return out
  end
  
 
# ---------------
# LOAD environment (libs and functions)
  if  occursin( r"size_structured", model_variation ) 
    fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )
  elseif  occursin( r"logistic_discrete", model_variation ) 
    fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  end
  

  include( fn_env )  # loads libs and setup workspace / data (fn_env is defined in the snowcrab_startup.jl)
     
