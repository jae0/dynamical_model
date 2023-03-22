
# starting environment for the DDE snow crab model


print( "WARNING: if this is the initial run, it will take a while to precompile/download all libraries" )


if false

  # if doing manual startup .. this is done automatically on start but in case it fails:

  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
  # include( joinpath( project_directory, "startup.jl" ))    # alt
 
  # translate model-specific functions, etc to generics
  # if model_variation=="logistic_discrete_basic"
  #   fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  # end

  # if model_variation=="logistic_discrete" 
  #   fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  # end

  # if model_variation=="logistic_discrete_map" 
  #   fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  
  # end

  fn_env = joinpath( project_directory, "logistic_discrete_environment.jl" )  

  include(  fn_env )
 

end


# load libs and check settings
# pkgs are defined in snowcrab_startup.jl
for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages


theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)

gr()

# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)

# allsavetimes = unique( vcat( survey_time, prediction_time  ) )

 
# fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_biomass.RData"
o = load( fndat, convert=true)
 
Y = o["Y"][∈(yrs).(o["Y"].yrs), :]
removals = o["L"][∈(yrs).(o["L"].yrs), :]

Kmu = [5.0, 60.0, 1.25]


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
        plot!( Y[:,:yrs] .+2, Y[:,:cfasouth_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfasouth_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfasouth_M4] )

        plot(  Y[:,:yrs], Y[:,:cfanorth_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfanorth_M1] )
        plot!( Y[:,:yrs] .+2 , Y[:,:cfanorth_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfanorth_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfanorth_M4] )


        plot(  Y[:,:yrs], Y[:,:cfa4x_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfa4x_M1] )
        plot!( Y[:,:yrs] .+2 , Y[:,:cfa4x_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfa4x_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfa4x_M4] )
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
  q = (  1.0, 0.1,  0.01, 10.0 ),
  qc = (-SminFraction, 0.1, -1.0, 1.0 ),
  m0 = ( 0.9, 0.2, 0.1, 1.0), 
  mlim =(0.0, 1.0),
  removed=removed,
  S = S,
  SminFraction = SminFraction,
  iok = iok,
  yeartransition = 0
) 


# translate model-specific functions, etc to generics
if model_variation=="logistic_discrete_historical"

  if (aulab=="cfanorth") | (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
  end

  PM = @set PM.K = (kmu, 0.25*kmu, kmu/10.0, kmu*10.0 )
  PM = @set PM.r = (1.0, 0.1, 0.25, 2.0)
  PM = @set PM.bpsd = (0, 0.1, 1.0e-9, 0.5) 
  PM = @set PM.bosd = (0, kmu*0.5, 1.0e-9, kmu/2.0)
  PM = @set PM.q = (1.0, 0.5,  1.0e-1, 10.0)
  PM = @set PM.m0 = ( 5, 5)
  PM = @set PM.mlim = ( 0.0, 1.25)
  
  fmod = logistic_discrete_turing_historical( PM )  # q only

elseif model_variation=="logistic_discrete_basic"
  
  if (aulab=="cfanorth") | (aulab=="cfasouth")
    PM = @set PM.yeartransition = 6
  end

  PM = @set PM.bpsd = ( 0.1, 0.05, 0.01, 0.25 )
  PM = @set PM.bosd = ( 0.1, 0.05, 0.01, 0.25 )
  PM = @set PM.q = (  1.0, 0.1,  0.5, 1.5 )
  
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
  PM = @set PM.q = (1.0, 0.1,  0.01, 10.0)
  PM = @set PM.qc = (0.0, 0.1, -1.0, 1.0)

  fmod = logistic_discrete_map_turing( PM )  

end
 

# run model estimations / overrides
Turing.setprogress!(false);
n_adapts=6000
n_samples=5000
n_chains=4


# NUTS-specific
# see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
rejection_rate = 0.65  ## too high and it become impossibly slow .. this is a good balance between variability and speed
max_depth = 7  ## too high and it become impossibly slow
init_ϵ = 0.01

if model_variation=="logistic_discrete_historical"   # pre-2022, mimic STAN defaults
  n_adapts=10000
  n_samples=5000
  rejection_rate = 0.99  
  max_depth=14  ## too high and it become impossibly slow
end

 
turing_sampler = Turing.NUTS(n_samples, rejection_rate; max_depth=max_depth, init_ϵ=init_ϵ )

print( model_variation, ": ", aulab, year_assessment )

res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )  

print( "results file:",  res_fn )


