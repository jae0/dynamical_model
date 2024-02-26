
# starting environment for the discrete logistic  snow crab model
 
if false

  # if doing manual startup .. this is done automatically on start but in case it fails:

  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
  
  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  

  include(  fn_env )
 
end


# # pkgs are defined in snowcrab_startup.jl and loaded in it
# for pk in pkgs
#   @eval using $(Symbol(pk))
# end   # Pkg.add( pkgs ) # add required packages

 

# plotting backend
# plotly(ticks=:native)                  # plotlyjs for richer saving options
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)
# gr(size = (300, 300), legend = false)  # provide optional defaults
gr()

fndat  = joinpath( bio_data_directory, "biodyn_biomass.RData" )
# fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_biomass.RData"
  
o = load( fndat, convert=true)
 
Y = o["Y"][∈(yrs).(o["Y"].yrs), :]
removals = o["L"][∈(yrs).(o["L"].yrs), :]

Kmu = [5.0, 65.0, 1.5]

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
        plot(  Y[:,:yrs], Y[:,:cfasouth] )
        plot(  Y[:,:yrs], Y[:,:cfasouth_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfasouth_M1] )
 
        plot(  Y[:,:yrs], Y[:,:cfanorth_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfanorth_M1] )
        
        
        plot(  Y[:,:yrs], Y[:,:cfa4x_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfa4x_M1] )
        
      end

 
nT = length(yrs)
nP = 5  # number of predictions into future (with no fishing)
nM = nP + nT  # total number of prediction years

dt = 1  # time resolution of solutions .. discrete annual so 1 year

nS = 1 # no state variables, not used 

no_digits = 3  # time floating point rounding 

smallnumber = 1.0e-9 # floating point value sufficient to assume 0 valued
    
# "survey index"
S = Y[:,Symbol("$aulab"  )]

# scale index where required
Smean = mean(skipmissing(S))
Sstd = std( skipmissing(S))
Smin = minimum(skipmissing(S))
Smax = maximum( skipmissing(S))
Srange = Smax - Smin 

SminFraction = Smin ./ Srange  # used as informative prior mean in some runs

if model_variation=="logistic_discrete_historical"
  # do nothing (no scaling)
elseif occursin( r"scale_center", model_variation ) 
  S = (S .- Smean ) ./ Sstd  # scale to std and center to 0 
else 
  S = (S .- Smin ) ./ Srange    # default is to scale (min,max) -> (0,1)
end


# id index
ki = aulab=="cfanorth" ? 1 :
     aulab=="cfasouth" ? 2 :
     aulab=="cfa4x"    ? 3 :
     0  # default

kmu = Kmu[ki]  

survey_time = Y[:,:yrs]   # time of observations for survey
# prediction_time = floor.(vcat( survey_time, collect(1:nP) .+ maximum(survey_time) ) )     
fish_time = round.( round.( removals[:,:ts] ./ dt; digits=0 ) .* dt; digits=no_digits)   # time of observations for landings

removed = removals[:,Symbol("$aulab")]
  
smallnumber = 1.0 / kmu / 10.0  # floating point value of sufficient to assume 0 valued

no_digits = 3  # time floating point rounding

survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)    # time of observations for survey
 
predtime =  4.0/12.0  # predictions ("m") are prefishery .. arbitrarily setting to 4/12
prediction_time =
  floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  #yrs
  round( predtime/dt ; digits=no_digits)   # april (m== prefishery)

iok = findall( !ismissing, S )
 
include( "logistic_discrete_functions.jl" )  #specific to model form


# basic params for "logistic_discrete"
PM = (
  yrs=yrs,
  nS = nS, 
  nT = length(yrs),
  nP = nP,  # number of predictions into future (with no fishing)
  nM = nM,  # total number of prediction years
  K = (kmu, 0.25*kmu, kmu/5.0, kmu*5.0 ),
  r = (1.0, 0.1, 0.5, 1.5),
  bpsd = ( 0.1, 0.05, 0.01, 0.5 ),
  bosd = ( 0.1, 0.05, 0.01, 0.5 ),
  q1 = (  1.0, 0.2,  0.01, 10.0 ),
  q0 = ( SminFraction, 0.1, -1.0, 1.0 ),
  m0 = ( 0.9, 0.2, 0.1, 1.0), 
  mlim =(0.0, 1.0),
  removed=removed,
  S = S,
  iok = iok,
  yeartransition = 0
) 



# run model estimations / overrides
Turing.setprogress!(false);

# Turing specific default options
n_adapts, n_samples, n_chains = 10000, 10000, 4
thin = 0

# Turing NUTS-specific default options  ..  see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
target_acceptance_rate, max_depth, init_ϵ = 0.65, 7, 0.01


# translate model-specific functions, etc to generics
if model_variation=="logistic_discrete_historical"

  if (aulab=="cfanorth")  
    PM = @set PM.yeartransition = 6
    PM = @set PM.K = kmu .* (1.0, 0.1, 1.0/4.0, 4.0 )
    PM = @set PM.r = (1.0, 0.1, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.bosd = kmu .* ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.q1 = ( 1.0, 0.5, 0.01, 2.0 )
    PM = @set PM.mlim =( 0.01, 1.25 )
    PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
     
    target_acceptance_rate, max_depth, init_ϵ = 0.65, 8, 0.01

  end

  if (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
    PM = @set PM.K = kmu .* (1.0, 0.1, 1.0/4.0, 4.0 )
    PM = @set PM.r = (1.0, 0.1, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.bosd = kmu .* ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.q1 = ( 1.0, 0.5, 0.01, 2.0 )
    PM = @set PM.mlim =( 0.01, 1.25 )
    PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
    
    target_acceptance_rate, max_depth, init_ϵ = 0.65, 8, 0.01

  end

  if (aulab=="cfa4x")  
    PM = @set PM.yeartransition = 0
    PM = @set PM.K = kmu .* (1.0, 0.1, 1.0/4.0, 4.0 )
    PM = @set PM.r = (1.0, 0.1, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.bosd = kmu .* ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.q1 = ( 1.0, 0.5, 0.01, 2.0 )
    PM = @set PM.mlim =( 0.01, 1.25 )
    PM = @set PM.m0 = (0.5, 0.25, 0.0, 1.25 )
     
    target_acceptance_rate, max_depth, init_ϵ = 0.65, 8, 0.01

  end
  
  
  fmod = logistic_discrete_turing_historical( PM )  # q1 only

elseif model_variation=="logistic_discrete_basic"
  # same as historical but normalize to reduce influce of large magnitudes and faster convergence

  if (aulab=="cfanorth")  
    PM = @set PM.yeartransition = 6
    PM = @set PM.K = (kmu, 0.25*kmu, kmu/4.0, kmu*4.0 )
    PM = @set PM.r = (1.0, 0.1, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.bosd = ( 0.1, 0.1, 0.01, 0.25 )
    PM = @set PM.q1 = (  1.0, 0.1,  0.25, 1.25 )
    PM = @set PM.mlim =( 0.1, 1.25 )
    PM = @set PM.m0 = ( 0.9, 0.2, 0.2, 1.25)

  end

  if (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
    PM = @set PM.K = (kmu, 0.25*kmu, kmu/4.0, kmu*4.0 )
    PM = @set PM.r = (1.0, 0.1, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.001, 0.25 )
    PM = @set PM.bosd = ( 0.1, 0.1, 0.001, 0.25 )
    PM = @set PM.q1 = (  1.0, 0.1,  0.25, 1.25 )
    PM = @set PM.mlim =( 0.2, 1.25 )
    PM = @set PM.m0 = ( 0.5, 0.2, 0.0, 1.25)

    target_acceptance_rate, max_depth, init_ϵ = 0.65, 7, 0.25

  end

  if (aulab=="cfa4x")  

    PM = @set PM.yeartransition = 0
    PM = @set PM.K = (kmu, 0.2*kmu, kmu/4.0, kmu*4.0 )
    PM = @set PM.r = (1.0, 0.2, 0.25, 1.75)
    PM = @set PM.bpsd = ( 0.1, 0.1, 0.001, 0.25 )
    PM = @set PM.bosd = ( 0.1, 0.1, 0.001, 0.25 )
    PM = @set PM.q1 = ( 1.0, 0.2, 0.25, 1.25 )
    PM = @set PM.mlim =( 0.1, 1.25 )
    PM = @set PM.m0 = ( 0.5, 0.2, 0.25, 1.25)

    target_acceptance_rate, max_depth, init_ϵ = 0.65, 7, 0.05

  end
  
  fmod = logistic_discrete_turing_basic( PM )   

elseif model_variation=="logistic_discrete"

  if (aulab=="cfanorth") | (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
  end

  fmod = logistic_discrete_turing( PM )   

elseif model_variation=="logistic_discrete_map"

  # not used .. just for testing
  if (aulab=="cfanorth") | (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
  end
  PM = @set PM.K = (kmu, 0.2*kmu, kmu/5.0, kmu*5.0 )
  PM = @set PM.r = (1.0, 0.1, 0.5, 3.0)
  PM = @set PM.bpsd = (0, 0.05, 0.01, 0.25) 
  PM = @set PM.bosd = (0, 0.05, 0.01, 0.25)
  PM = @set PM.q1 = (1.0, 0.1,  0.01, 10.0)
  PM = @set PM.q0 = (0.0, 0.1, -1.0, 1.0)

  fmod = logistic_discrete_map_turing( PM )  

end
 

if model_variation=="logistic_discrete_historical"   # pre-2022, mimic STAN defaults
  n_adapts, n_samples, n_chains = 10000, 10000, 4
  target_acceptance_rate =  0.65  #0.99
  max_depth=  7  # 14  ## too high and it become impossibly slow
  thin=0
end


# by default use NUTS sampler ... SMC is another good option if NUTS is too slow
turing_sampler = Turing.NUTS(n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )
print( model_variation, ": ", aulab, year_assessment )

# mcmc save file name and location
res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )  
print( "results file:",  res_fn )


