
# starting environment for the DDE snow crab model


if false

  # if doing manual startup .. this is done automatically on start but in case it fails:

  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
   
  fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )  
  include(  fn_env )
    
end


# load libs and check settings


# add Turing@v0.21.10  # to add a particular version

pkgs = [
  "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random",
  # "DynamicHMC", "AdvancedHMC",  "AdvancedMH",  "DynamicPPL",  "AbstractPPL",  "Memoization",
  "ForwardDiff", "DataFrames", "JLD2", "PlotThemes", "Colors", "ColorSchemes", "RData",
  "Plots", "StatsPlots",  "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
  "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra"
]

for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages


theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)

gr()

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
 

# perpare dat for dde run of fishery model

fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
o = load( fndat, convert=true)
Y = o["Y"]

Kmu = o["Kmu"]

Kmu = [5.5, 60.0, 1.5]

removals = o["L"]
MW = o["M0_W"]


    if false
        # alternatively, if running manually:
        # can run R-code that creates local RData file with required data
        # run in R externally or from within julia or ..
        # from within julia

        # using RCall
        # # type $ in Julia's command prompt starts an R session.
        # # .. run below
        # # type <backspace> to escape back to julia

        # source( file.path( code_root, "bio_startup.R" )  )
        # require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        # fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/12)

        # # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
        # @rget Y
        # @rget Kmu
        # @rget removals
        # @rget ty

        # example line plots
        plot(  Y[:,:yrs], Y[:,:cfasouth_M0] )
        plot!( Y[:,:yrs]  , Y[:,:cfasouth_M1] )
        plot!( Y[:,:yrs] , Y[:,:cfasouth_M2] )
        plot!( Y[:,:yrs] , Y[:,:cfasouth_M3] )
        plot!( Y[:,:yrs] , Y[:,:cfasouth_M4] )

        plot(  Y[:,:yrs], Y[:,:cfanorth_M0] )
        plot!( Y[:,:yrs]  , Y[:,:cfanorth_M1] )
        plot!( Y[:,:yrs]  , Y[:,:cfanorth_M2] )
        plot!( Y[:,:yrs] , Y[:,:cfanorth_M3] )
        plot!( Y[:,:yrs] , Y[:,:cfanorth_M4] )


        plot(  Y[:,:yrs], Y[:,:cfa4x_M0] )
        plot!( Y[:,:yrs]  , Y[:,:cfa4x_M1] )
        plot!( Y[:,:yrs]  , Y[:,:cfa4x_M2] )
        plot!( Y[:,:yrs] , Y[:,:cfa4x_M3] )
        plot!( Y[:,:yrs] , Y[:,:cfa4x_M4] )
    end


nT = length(yrs)
nP = 5  # number of predictions into future (with no fishing)
nM = nP + nT  # total number of prediction years

nS = 6  # no. state variables


# "survey index"
statevars = [
  Symbol("$aulab","_M0"),
  Symbol("$aulab","_M1"),
  Symbol("$aulab","_M2"),
  Symbol("$aulab","_M3"),
  Symbol("$aulab","_M4"),
  Symbol("$aulab","_f_mat")
]

S = Matrix(Y[:, statevars ])
Smins = [minimum(skipmissing(S[:,i])) for i in 1:nS ]
Smaxs = [maximum(skipmissing(S[:,i])) for i in 1:nS ]

# scale index to min-max
for i in 1:nS
  S[:,i] = (S[:,i] .- Smins[i] ) ./ ( Smaxs[i] .- Smins[i] ) # range from 0=min to 1=max
end

# interpolating function for mean weight
mwspline = extrapolate( interpolate( MW[:,Symbol("mw_", "$aulab") ], (BSpline(Linear()) ) ),  Interpolations.Flat() )
mw = Interpolations.scale(mwspline, yrs )

scale_factor = mw(yrs) / (1000 *1000 ) # convert numbers to kt biomass , also used in plots


# convert to (biomass kt to number)

# id index
ki = aulab == "cfanorth" ? 1 :
     aulab == "cfasouth" ? 2 :
     aulab == "cfa4x"    ? 3 :
     0  # default

kmu  =  Kmu[ki] / mean(scale_factor)

smallnumber = 1.0 / kmu / 10.0  # floating point value of sufficient to assume 0 valued

no_digits = 3  # time floating point rounding


dt =  aulab == "cfanorth" ? 0.05 :
      aulab == "cfasouth" ? 0.05 :
      aulab == "cfa4x"    ? 0.05  :   # long fishing seasons .. aggregate landings more stable
      0.1   # default

# spin up time of ~ 1 cycle prior to start of dymamics and project nP years into the future
tspan = (minimum(yrs) - 8.1, maximum(yrs) + nP + 1.1 )


survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)    # time of observations for survey
Si = findall( x-> !ismissing(x), vec(sum(S, dims=2)))  # compute data likelihoods only when data exist ... to speed up comps
nSI = length(Si)

# this only adds habitat space  ... predation is also a useful one ..
# speed is the issue
predtime =  9.0/12.0
prediction_time =
  floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  #yrs
  round(round( predtime/dt; digits=0 ) *dt; digits=no_digits)   # sept

#  sa to fraction
external_forcing =  reshape( [
   Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
   Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
   Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
   Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
   Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
   Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] )
  ], nT, nS )


efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
hsa = Interpolations.scale(efc, yrs .+ predtime, 1:nS )

fish_time =  round.( round.( removals[:,:ts]  ./ dt; digits=0 ) .* dt; digits=no_digits)    # time of observations for survey

ys =  aulab == "cfanorth" ? "yrs" :
      aulab == "cfasouth" ? "yrs" :
      aulab == "cfa4x"    ? "yrs_4x"  :
      "yrs"   # default

fish_year =  round.( round.( removals[:,Symbol(ys)] ./ dt; digits=0 ) .* dt; digits=no_digits)    # time of observations for survey
 
removed = removals[:,Symbol("$aulab")]

# keep nonzero elements
ikeep = findall( x-> x>0, removed)
if length(ikeep) > 0
  removed = removed[ ikeep ]
  fish_time = fish_time[ikeep]
  fish_year = fish_year[ikeep]
end
  
directory_output = joinpath( project_directory, "outputs", model_variation )
mkpath(directory_output)

# model-specifics functions and data

include( "fishery_model_functions.jl" )  # to load core dynamical model functions


if model_variation=="size_structured_dde_normalized" 
  include( "size_structured_dde_normalized_functions.jl" )  

  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ]  ; # generics to bootstrap the process

  # history function (prior to start)  defaults to values of 0.5*u before t0; note u = (0,1)
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS)  .* 0.9

  function affect_fishing!(integrator)
    i = findall(t -> t == integrator.t, fish_time)[1]
    integrator.u[1] -=  removed[ i ] / integrator.p[2][1]  # p[2] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
  end
   
elseif  model_variation=="size_structured_dde_unnormalized"
  include( "size_structured_dde_unnormalized_functions.jl" )  

  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ]  .* kmu; # generics to bootstrap the process

  # history function (prior to start)  defaults to values of 0.5*u before t0; note u = (0,1)
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS)  .* 0.9 .* kmu

  function affect_fishing!(integrator)
    i = findall(t -> t == integrator.t, fish_time)[1]
    integrator.u[1] -=  removed[ i[1] ]  # sol on same scale
  end
   
   
else 
  
  error("model_variation not found")

end
 

# callbacks for external perturbations to the system (deterministic fishing without error)
cb = PresetTimeCallback( fish_time, affect_fishing! )
# alternative formulation:
# cb = CallbackSet(
#   PresetTimeCallback( fish_time, affect_fishing! ),
#   PositiveDomain()
# );


# delay resolution
tau = [1.0]  

# Turing sampling-specific
n_adapts=1000
n_samples=1000
n_chains=4

# NUTS-specific
rejection_rate = 0.75
max_depth=10
init_ϵ=0.01

# DiffEq-specific
p = dde_parameters() # dummy values needed to bootstrap DifferentialEquations/Turing initialization
prob = DDEProblem( size_structured_dde!, u0, h, tspan, p, constant_lags=tau  )  # tau=[1]

# Turing-DiffEq model
fmod = size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM, solver, dt )

