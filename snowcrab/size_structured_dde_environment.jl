
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


theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)


# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)
gr()

# prepare dat for dde run of fishery model
print( "\n\nThe following RData warnings can be safely ignored. \n\n" )

fndat  = joinpath( bio_data_directory, "biodyn_number_size_struct.RData" )
 
o = load( fndat, convert=true)

Yyrs = floor.(Int, o["Y"].yrs)
Y = o["Y"][∈(yrs).(Yyrs), :]

removalsyrs = floor.(Int, o["L"].yrs)
removals = o["L"][∈(yrs).(removalsyrs), :]  # in numbers (not mass)

MW = o["M0_W"][∈(yrs).(o["M0_W"].mw_yrs), :]
MW.yrs = MW.mw_yrs


Kmu = [5.5, 60.0, 1.25]   # 5.0, 60.0, 1.25 

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
 
        plot(  Y[:,:yrs], Y[:,:cfanorth_M0] )
        plot!( Y[:,:yrs]  , Y[:,:cfanorth_M1] )
  

        plot(  Y[:,:yrs], Y[:,:cfa4x_M0] )
        plot!( Y[:,:yrs]  , Y[:,:cfa4x_M1] )
  
    =#


nT = length(yrs)
nP = 5  # number of predictions into future (with no fishing)
nM = nP + nT  # total number of prediction years

nS = 6  # no. state variables

# BLY =[8]
BLY = [7,8,9] #  birth lag years, centered on 8 yrs offset

nB = length(BLY)

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

SminFraction = Smin ./ Srange  # used as informative prior mean in some runs

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
     
dt = (0.02, 0.02, 0.02)[ki]    # resolution of time (fraction of year)  (i.e. days = dt*365; weeks=dt*52); dt=0.02 is ~ weekly

# spin up time of ~ 1 cycle prior to start of dymamics and project nP years into the future
tspan = (minimum(yrs) - 11.1, maximum(yrs) + nP + 1.1 )

survey_time = discretize_decimal( Y[:,:yrs], dt)     

Si = findall( x-> !ismissing(x), vec(sum(S, dims=2)))  # compute data likelihoods only when data exist ... to speed up comps
nSI = length(Si)

# this only adds habitat space  ... predation is also a useful one ..
predtime =  discretize_decimal( 9.0/12.0, dt )
prediction_time = floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+ predtime   # sept

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

ys = ( "yrs", "yrs", "yrs_4x")[ki]

# fishing pattern seasonal over past 5 yrs (function) .. used for projections
fish_time_max = discretize_decimal( maximum(removals[:,:ts]), dt)  + dt 
fp_time = fish_time_max:dt:maximum(prediction_time)
fish_time_project = discretize_decimal( collect(fp_time), dt )     

# fish_year is fishery "year" 
fsh = DataFrame(
  :fish_time => discretize_decimal( removals[:,:ts] , dt) ,
  :fish_year => discretize_decimal( removals[:,Symbol(ys)], dt),
  :removed => removals[:,Symbol("$aulab")]  #number
)

# keep nonzero elements
fsh = fsh[findall( x-> x>0, fsh.removed),:]

# aggregate if required
afsh = DataFrame( combine( groupby(fsh, [:fish_time]), :removed => sum, :fish_year => unique) )

# make available as a global variable
fish_time = afsh.fish_time
fish_year = afsh.fish_year_unique
removed = afsh.removed_sum
 
# model-specifics functions and data
 
# DiffEq-model setup

if model_variation=="size_structured_dde_normalized" 
  include( "size_structured_dde_normalized_functions.jl" )  
  function affect_fishing!(integrator)
    i = findall(t -> t == integrator.t, fish_time)[1]
    integrator.u[1] -=  removed[ i ] / integrator.p[3][1]  # p[3] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
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
# cb = CallbackSet( PresetTimeCallback( fish_time, affect_fishing! ), PositiveDomain() );
 
# generics to bootstrap the process; initial conditions .. "pristine", unfished condition ... 
# however, natural mortality might have been different so be careful here
u0 = 
  model_variation=="size_structured_dde_unnormalized" ? ones(nS) .* kmu .* 1.0 :
  model_variation=="size_structured_dde_normalized"   ? ones(nS) .* 1.0 :
  1.0

# history function (prior to start)  defaults to values of kmu / 1 before t0;  
# h(p, t; idxs=nothing) = typeof(idxs) <: Number ? u0 : u0 
h(p, t) = u0 

tau = BLY  # delay resolution
p = dde_parameters() # dummy values needed to bootstrap DifferentialEquations/Turing initialization
prob = DDEProblem{true}( size_structured_dde!, u0, h, tspan, p, constant_lags=tau  )  #  , neutral=true create container for problem definition 

# choose DiffEq solver:
# stiff solvers: Rodas4()  ; Rosenbrock23()
# diffeq_solver = MethodOfSteps(Rosenbrock23()) # slow
diffeq_solver = MethodOfSteps(AutoTsit5(Rodas5()))
# diffeq_solver = MethodOfSteps(AutoTsit5(KenCarp47()))
# diffeq_solver = MethodOfSteps(KenCarp47())  
# diffeq_solver = MethodOfSteps(Tsit5())   # faster
# diffeq_solver = MethodOfSteps( Rodas4() )
# diffeq_solver = MethodOfSteps(Rodas5())  # safer
 
solver_params = (
  prob=prob,
  abstol = 1.0e-9,    
  reltol = 1.0e-9, 
  dt = dt,
  saveat = collect(tspan[1]:dt:tspan[2]),
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

# pl = plot(x->pdf(LogNormal(log(10), 1.0), x), xlim=(0,10)) #

# choose model and over-rides if any
 
PM = (
    nS = nS, 
    nSI = nSI,
    nB = nB,
    nG = 4,  # n transition moults .. growth
    nT = length(yrs),
    nP = 5,  # number of predictions into future (with no fishing)
    nM = nP + nT,  # total number of prediction years
    logkmu = (logkmu, 0.1),
    logScv = (logScv, 0.1),
    b = ( log(1.0), 0.1),
    d =  ( log( exp(0.25)-1.0 ), 0.1 ),
    d2 = ( log( exp(0.75)-1.0 ), 0.1 ),
    v =  ( log( exp(0.90)-1.0 ), 0.1 ),
    q1 = ( 0.5, 0.1 ),  # assume about 50% overly optimistic estimation from CARSTM (due to assumption that each areal unit is homogenous)
    q0 = ( SminFraction ./ 2.0, 0.1 ),  # lower detection limit of survey/CARSM
    BLY = BLY,
    Si = Si,
    S = S,
    data = S[Si,:],
    # datavector = vec(S[Si,:]),  # for MVN .. no advantage
    Stime = survey_time[Si],
    eps=1e-9
) 
# the following tweak Lognormal priors by area  


# Turing specific default options
n_adapts, n_samples, n_chains = 3000, 5000, 4

# Turing NUTS-specific default options  ..  see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
target_acceptance_rate, max_depth, init_ϵ = 0.65, 7, 0.01

# model-specific overrides
if model_variation=="size_structured_dde_normalized" 
   
    if aulab=="cfanorth"

        PM = @set PM.b =  ( log(2.0), 0.1 ) 
        
        PM = @set PM.q0 = ( SminFraction ./ 2.0, 0.1)
        PM = @set PM.q1 = ( 0.5, 0.1 )
  
    elseif aulab=="cfasouth" 
        
        # higher birth rates expected due to temperaturess
        PM = @set PM.b =  ( log(2.0), 0.25 )

        PM = @set PM.q0 = ( SminFraction ./ 2.0, 0.15)
        PM = @set PM.q1 = ( 0.5, 0.15 )
        
    elseif aulab=="cfa4x" 

        # higher birth rates expected due to temperaturess
        PM = @set PM.b =  ( log(2.0), 0.1 )

        # suggests CARSTM is less biased in 4X
        PM = @set PM.q0 =  ( SminFraction, 0.1 )
        PM = @set PM.q1 =  ( 1.0, 0.1 )
  
    end

    # define model
    fmod = size_structured_dde_turing( PM, solver_params )

elseif  model_variation=="size_structured_dde_unnormalized"
     
    print( "warning: model needs some updating .. do not use until it is checked" )
    
    fmod = size_structured_dde_turing( S, kmu, tspan, prob, nS )

end
 
 
# by default use SMC sampler ... SMC is another good option if NUTS is too slow
# turing_sampler = Turing.NUTS(n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )
turing_sampler = Turing.SMC()
print( string( model_variation, " : ", aulab, " - ", year_assessment) )

# mcmc save file name and location
res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )  
print( "results file:",  res_fn )

