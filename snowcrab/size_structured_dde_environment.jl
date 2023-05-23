
# starting environment for the DDE snow crab model


print( "WARNING: if this is the initial run, it will take a while to precompile/download all libraries" )


if false

  # if doing manual startup .. this is done automatically on start but in case it fails:

  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
   
  fn_env = joinpath( project_directory, "size_structured_dde_environment.jl" )  
  include(  fn_env )
    
end




theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)

gr()

# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)

# allsavetimes = unique( vcat( survey_time, prediction_time  ) )

# prepare dat for dde run of fishery model
# ---------------
# make a copy of the input data in case ... 

fndat_source = joinpath( bio_data_directory, "bio.snowcrab", "modelled", 
    "1999_present_fb", "fishery_model_results", "turing1", "biodyn_number_size_struct.RData" )

fndat = joinpath( model_outdir, basename(fndat_source) )

if (!isfile(fndat)) 
  # prompt to input
  print("\nData file not found. Copy from: \n")
  print(fndat_source)
  print("\nTo: \n")
  print( fndat )
  print( "\nType 'Yes' to proceed >  ")
  confirm = readline()
  if confirm=="Yes"
    cp( fndat_source, fndat; force=true )
  end
end

o = load( fndat, convert=true)

Yyrs = floor.(Int, o["Y"].yrs)
Y = o["Y"][∈(yrs).(Yyrs), :]

removalsyrs = floor.(Int, o["L"].yrs)
removals = o["L"][∈(yrs).(removalsyrs), :]  # in numbers (not mass)

MW = o["M0_W"][∈(yrs).(o["M0_W"].mw_yrs), :]
MW.yrs = MW.mw_yrs


Kmu = [5.0, 65.0, 1.25]   # 5.0, 60.0, 1.25 

    #=
        # alternatively, if running manually:
        # can run R-code that creates local RData file with required data
        # run in R externally or from within julia or ..

        # from within julia

        using RCall
        # type $ in Julia's command prompt starts an R session.
        # .. run below
        # type <backspace> to escape back to julia

        source( file.path( code_root, "bio_startup.R" )  )
        require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        fishery_model_data_inputs( year.assessment=year.assessment, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/12)

        # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
        @rget Y
        @rget Kmu
        @rget removals
        @rget ty

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
    =#


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


# scale index where required
Smean = [mean(skipmissing(S)) for i in 1:nS ]
Sstd = [std( skipmissing(S)) for i in 1:nS ]
Smin = [minimum(skipmissing(S[:,i])) for i in 1:nS ]
Smax = [maximum(skipmissing(S[:,i])) for i in 1:nS ]
Srange = Smax .- Smin 

SminFraction = - Smin ./ Srange  # used as informative prior mean in some runs

# CV of Y
statevars_sd = [
  Symbol("$aulab", "_sd", "_M0"),
  Symbol("$aulab", "_sd", "_M1"),
  Symbol("$aulab", "_sd", "_M2"),
  Symbol("$aulab", "_sd", "_M3"),
  Symbol("$aulab", "_sd", "_M4"),
  Symbol("$aulab", "_sd", "_f_mat")
]

Ssd = Matrix(Y[:, statevars_sd ])

Scv = Ssd ./ S  

# scale index to min-max
if occursin( r"unnormalized", model_variation )
  # do nothing (no scaling)
elseif occursin( r"scale_center", model_variation ) 
  for i in 1:nS
    S[:,i] = (S[:,i] .- Smean[i] ) ./ Sstd[i]    # scale to std and center to 0 
  end
else 
  # default is to normalize (min, max) to (0,1)
  for i in 1:nS
    S[:,i] = (S[:,i] .- Smin[i] ) ./ Srange[i]   # range from 0=min to 1=max
  end
end

# deal with missing CV's
for i in 1:nS
  u = findall( x-> ismissing(x), Scv[:,i] )  
  if length(u) > 0
    Scv[u,i] .= S[u,i] # ie. poisson
  end
end
for i in 1:nS
  u = findall( x-> ismissing(x), Scv[:,i] )  
  if length(u) > 0
    Scv[u,i] .= 0.5 # ie. no information ( u is in interval 0,1 .. 0.5 covers it nicely )
  end
end

logScv = log.(Scv) # on log scale to reduce further computations


# interpolating function for mean weight (of fb only)
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
logkmu = log(kmu)

smallnumber = 1.0 / (kmu * 10.0) # floating point value of sufficient to assume 0 valued
     
no_digits = 3  # time floating point rounding

dt = (0.01, 0.01, 0.01)[ki] 

# spin up time of ~ 1 cycle prior to start of dymamics and project nP years into the future
tspan = (minimum(yrs) - 10.1, maximum(yrs) + nP + 1.1 )


survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)     
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

fish_time =  round.( round.( removals[:,:ts]  ./ dt; digits=0 ) .* dt; digits=no_digits)     

ys = ( "yrs", "yrs", "yrs_4x")[ki]

fish_year =  round.( round.( removals[:,Symbol(ys)] ./ dt; digits=0 ) .* dt; digits=no_digits)   # fishery "year"  
 
removed = removals[:,Symbol("$aulab")]  #number

# keep nonzero elements
ikeep = findall( x-> x>0, removed)
if length(ikeep) > 0
  removed = removed[ ikeep ]
  fish_time = fish_time[ikeep]
  fish_year = fish_year[ikeep]
end
  
# fishing pattern seasonal over past 5 yrs (function) .. used for projections
fish_time_max = maximum(removals[:,:ts]) + dt 
fp_time = fish_time_max:dt:maximum(prediction_time)
fish_time_project =  round.( round.( fp_time ./ dt; digits=0 ) .* dt; digits=no_digits)     



# model-specifics functions and data


# DiffEq-model setup

if model_variation=="size_structured_dde_normalized" 
  include( "size_structured_dde_normalized_functions.jl" )  
 
  function affect_fishing!(integrator)
    i = findall(t -> t == integrator.t, fish_time)[1]
    integrator.u[1] -=  removed[ i ] / integrator.p[2][1]  # p[2] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
  end
 

elseif  model_variation=="size_structured_dde_unnormalized"
  include( "size_structured_dde_unnormalized_functions.jl" )  
 
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



# generics to bootstrap the process; initial conditions
u0 = 
  model_variation=="size_structured_dde_unnormalized" ? ones(nS) .* kmu .* 0.5 :
  model_variation=="size_structured_dde_normalized"   ? ones(nS) .* 0.5 :
  0.5

# history function (prior to start)  defaults to values of kmu / 1 before t0;  
# h(p, t; idxs=nothing) = typeof(idxs) <: Number ? u0 : u0 
h(p, t) = u0 

tau = [1.0]  # delay resolution
p = dde_parameters() # dummy values needed to bootstrap DifferentialEquations/Turing initialization
prob = DDEProblem{true}( size_structured_dde!, u0, h, tspan, p, constant_lags=tau  )  # create container for problem definition 

#= 
    # choose DiffEq solver:
    # stiff solvers: Rodas4()  ; Rosenbrock23()
    diffeq_solver = MethodOfSteps(Rosenbrock23()) # slow
    diffeq_solver = MethodOfSteps(Rodas4())
    
    # other solvers: BS3() and Vern6() also RK4()
    diffeq_solver = MethodOfSteps(BS3())   
    diffeq_solver = MethodOfSteps(Rodas4()) 
    diffeq_solver = MethodOfSteps(Rosenbrock23()) 
    diffeq_solver = MethodOfSteps(Vern6())  
    diffeq_solver = MethodOfSteps(RK4())   
    diffeq_solver = MethodOfSteps(TRBDF2())  
    diffeq_solver = MethodOfSteps(QNDF())  
    diffeq_solver = MethodOfSteps(Vern7())  
    diffeq_solver = MethodOfSteps(KenCarp4())  
    
    diffeq_solver = MethodOfSteps(Tsit5())   # faster
    diffeq_solver = MethodOfSteps(Rodas5())  # safer
=#  

diffeq_solver = MethodOfSteps(AutoTsit5(Rosenbrock23()))

solver_params = (
  prob=prob,
  abstol = 1.0e-6,  # these are diffeq defaults 
  reltol = 1.0e-6, 
  dt = dt,
  # saveat = collect(tspan[1]:dt:tspan[2]),
  h = h,
  cb = cb,
  tspan = tspan,
  solver = diffeq_solver
)

#=  lognormal rates:
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    log( exp(0.5)-1.0 ) = log(0.6487 ) = -0.4328 .. 50%
  
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.6912 = 0.5254   .. ~99% (moult) transition rate per year (mode) 
=#

# choose model and over-rides if any
PM = (
    nS = nS, 
    nSI = nSI,
    nB = 2,
    nG = 4,  # n transition moults .. growth
    nT = length(yrs),
    nP = 5,  # number of predictions into future (with no fishing)
    nM = nP + nT,  # total number of prediction years
    logkmu = (logkmu, 0.25),
    logScv = (logScv, 0.25),
    b = ( log(1), 0.25),
    d =  ( log( exp(0.2)-1.0 ), 0.25 ),
    d2 = ( log( exp(0.5)-1.0 ), 0.5 ),
    v =  ( log( exp(0.9)-1.0), 0.5 ),
    q = (1.0, 0.1),
    qc = (SminFraction, 0.1),
    Si = Si,
    S = S,
    data = S[Si,:],
    # datavector = vec(S[Si,:]),  # for MVN .. no advantage
    Stime = survey_time[Si]
) 
# the following tweak Lognormal priors by area  

if model_variation=="size_structured_dde_normalized" 
    n_adapts=250
    n_samples=500
    n_chains=4
  
    rejection_rate = 0.65
    max_depth = 7
    init_ϵ = 0.05

    if aulab=="cfanorth"

        # PM = @set PM.b =  ( log(10), 0.50 ) 
        # PM = @set PM.d =  ( log( exp(0.2)-1.0 ), 0.25 ) 
        # PM = @set PM.d2 = ( log( exp(0.5)-1.0 ), 0.50 )
        # PM = @set PM.v =  ( log( exp(0.9)-1.0 ), 0.50 )

    elseif aulab=="cfasouth" 
      
        # solver_params = @set solver_params.abstol = 1.0e-9
        # solver_params = @set solver_params.reltol = 1.0e-9
        
        # PM = @set PM.b =  ( log(10), 0.50 ) 
        # PM = @set PM.d =  ( log( exp(0.2)-1.0 ), 0.25 ) 
        PM = @set PM.d2 = ( log( exp(0.4)-1.0 ), 0.50 )
        # PM = @set PM.v =  ( log( exp(0.8)-1.0 ), 0.50 )   #         testing

    elseif aulab=="cfa4x" 

        # PM = @set PM.b =  ( log(10), 0.50 ) 
        # PM = @set PM.d =  ( log( exp(0.2)-1.0 ), 0.25 ) 
        # PM = @set PM.d2 = ( log( exp(0.5)-1.0 ), 0.50 )
        # PM = @set PM.v =  ( log( exp(0.95)-1.0 ), 0.50 )

    end

    # define model
    fmod = size_structured_dde_turing( PM=PM, solver_params=solver_params )

elseif  model_variation=="size_structured_dde_unnormalized"

    n_adapts=1000
    n_samples=1000
    n_chains=4
    
    
    # NUTS-specific run options
    # see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
    rejection_rate = 0.65  ## too high and it become impossibly slow .. this is a good balance between variability and speed
    max_depth=7  ## too high and it become impossibly slow
    init_ϵ=0.0125
    
    print( "warning: model needs some updating" )
    
    # turing_sampler = Turing.NUTS(n_samples, rejection_rate ) #; max_depth=max_depth, init_ϵ=init_ϵ )
    fmod = size_structured_dde_turing( S, kmu, tspan, prob, nS )

end
 

# Turing sampling-specific run options (NUTS-specific)
# see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz

if model_variation=="size_structured_dde_normalized" 
    n_adapts=500
    n_samples=500
    n_chains=4
    rejection_rate = 0.65
    max_depth = 7
    init_ϵ = 0.01
elseif  model_variation=="size_structured_dde_unnormalized"
    n_adapts=1000
    n_samples=1000
    n_chains=4
    rejection_rate = 0.65  ## too high and it become impossibly slow .. this is a good balance between variability and speed
    max_depth=7  ## too high and it become impossibly slow
    init_ϵ=0.01
end

turing_sampler = Turing.NUTS(n_samples, rejection_rate; max_depth=max_depth, init_ϵ=init_ϵ )

print( string( model_variation, " : ", aulab, " - ", year_assessment) )


res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )  

print( "results file:",  res_fn )

