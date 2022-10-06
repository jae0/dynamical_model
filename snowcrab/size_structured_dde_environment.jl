
# starting environment for the DDE snow crab model

 
# load libs and check settings

# add Turing@v0.21.10  # to add a particular version
 
pkgs = [ 
  "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random",   
  # "DynamicHMC", "AdvancedHMC",  "AdvancedMH",  "DynamicPPL",  "AbstractPPL",  "Memoization", 
  "ForwardDiff", "DataFrames", "FillArrays", "JLD2", "PlotThemes", "Colors", "ColorSchemes",
  "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
  "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra" 
]

for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages
    

theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)

# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)
  
# allsavetimes = unique( vcat( survey_time, prediction_time  ) )

# stiff solvers: Rodas4()  ; Rosenbrock23()
# solver = MethodOfSteps(Rosenbrock23()) # slow  
# solver = MethodOfSteps(Rodas4())  
# other solvers: BS3() and Vern6() also RK4()
# solver = MethodOfSteps(Rodas5())  # safer 

# relative timings:
# solver = MethodOfSteps(Tsit5())  # 10 - 41.43 
# solver = MethodOfSteps(Rodas5())   # 20.94  - 71.73
# solver = MethodOfSteps(BS3())   # 56.1
# solver = MethodOfSteps(Rodas4()) #   24.86- 82.79
# solver = MethodOfSteps(Rosenbrock23()) #  71.48
# solver = MethodOfSteps(Vern6())  # 73.98
# solver = MethodOfSteps(RK4())   # 76.28
# solver = MethodOfSteps(TRBDF2())  # 92.16
# solver = MethodOfSteps(QNDF())  # 110.79
# solver = MethodOfSteps(Vern7())  #  111.7
# solver = MethodOfSteps(KenCarp4())  # 139.88


solver = MethodOfSteps(Tsit5())   # faster
# solver = MethodOfSteps(Rodas5())  # safer 

Turing.setadbackend(:forwarddiff)  # only AD that works right now

# Turing.setrdcache(true) # reverse diff not working right Newton
 

if false

    # if doing manual startup .. this is done automatically on start but in case it fails:
    
    project_directory = @__DIR__() #  same folder as the current file
    push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
    include( "startup.jl" )    
    # include( joinpath( project_directory, "startup.jl" ))    # alt
 
end

