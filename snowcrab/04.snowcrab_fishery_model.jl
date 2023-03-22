
# ------------------------------
# fishery model

# NOTE: this is JULIA (https://julialang.org/) code and not R 

# it requires a number of libraries for itself and this environment is set up in associated *.environment.jl files

# for more info on the main modelling tools:
# https://turing.ml/stable/
# https://docs.sciml.ai/Overview/stable/
# https://github.com/SciML/DifferentialEquations.jl


# DATA REQUIREMENTS: 
# if not already created in their respective scripts, the following R-code can be run (in R):

    #= 
    # (NOTEL do not alter the "#=" and "=#", these are julia block comment identifiers)
    # R-code ==== start
        
      # NOTE::: this requires 03.snowcrab_carstm.r to be completed
      source( file.path( code_root, "bio_startup.R" )  )
      loadfunctions("bio.snowcrab")
      
      # prep data for discrete version
      if (grepl("logistic_discrete", model_variation )) {
          fishery_model_data_inputs( year.assessment=year.assessment, type="biomass_dynamics", for_julia=TRUE   )
      }

      # prep data for continuous version: 
      if (grepl("size_structured", model_variation)) {
          # fishery landings has a weekly time step = 2/52 ~ 0.0385 ~ 0.04  X dt=0.01 seems to work best
          fishery_model_data_inputs( year.assessment=year.assessment, type="size_structured_numerical_dynamics", for_julia=TRUE, time_resolution=2/52  )
      }

    # R-code ==== end
    =#  

# DEFINE KEY DIRECTORIES:
  # these need to be defined  ... if not in start up call or ".julia/config/startup.jl" or  local startup.jl
  project_directory = joinpath( homedir(), "bio", "bio.snowcrab", "inst", "julia" ) 
  bio_data_directory = joinpath( homedir(), "bio.data" )  
  outputs_directory = joinpath( homedir(), "bio.data", "bio.snowcrab", "fishery_model" ) 
  

# RUN LEVEL OPTIONS: time, area model 
  year_assessment = 2022   # <<<<<<<<-- change
  yrs = 1999:year_assessment  
  
  model_variation = "size_structured_dde_normalized"  #      24hrs, >24 hrs

    #= 
        model_variation = "logistic_discrete_historical"    #   pre-2022 method  :: ~ 1 hr # pre-2022, no normalization, q-based observation model
        model_variation = "logistic_discrete_basic"     # q catchability only for observation model
        model_variation = "logistic_discrete"              # q and intercept for observation model
        model_variation = "logistic_discrete_map"  # logistic map ... more extreme fluctuations
        model_variation = "size_structured_dde_normalized"  #      24hrs, >24 hrs
        model_variation = "size_structured_dde_unnormalized"  # basic continuous model without normaliztion ... very very slow .. do not use
    (incomplete params need tweaking) ::   24hrs, >24 hrs
    =#
  
  # choose a region of interest"
  aulab ="cfanorth"   
  aulab ="cfasouth"   
  aulab ="cfa4x"     

  
# ---------------
# define global variables and model-specific save location
  include( joinpath(project_directory, "snowcrab_startup.jl" ) )  # add some paths and package requirements

    #= 
        # should this be a first time run, to install libs:
        # this is done automatically on first run (if file is not copied)
        # cp( fndat_source, fndat; force=true )   # to force copy of data file 
        for pk in pkgs; Pkg.add(string(Symbol(pk))); end   
        
        # or a single package at a time:, etc
        Pkg.add( "Turing" ) 
    =#
    


# ---------------
# Development & testing  -- can safely skip over

#   determine params that have reasonable distribution and extract modes/means
#   repeat (manually) until found (record random number seed to have direct control )         
#   loads libs and setup workspace / data (fn_env is defined in the snowcrab_startup.jl)
#   include( joinpath(project_directory, "snowcrab_startup.jl" ) )  # add some paths and package requirements


    #=
        # debugging/development: to test dynamical model with generic/random parameters
        (test, pl) = fishery_model_test( "basic" ); pl
        (test, pl) = fishery_model_test( "random_external_forcing"  ); pl
        (test, pl) = fishery_model_test( "fishing"  ); pl
        (test, pl) = fishery_model_test( "nofishing" ); pl
        showall( summarize( test ) )
      
        
        # CHoose AD backend -- forwarddiff is stable
        using SciMLSensitivity
        using ForwardDiff
        using Turing
        Turing.setadbackend(:forwarddiff)   
        using ReverseDiff # slow
        Turing.setadbackend(:reversediff)  # only AD that works right now
        Turing.setrdcache(true)
        using Zygote # ok
        Turing.setadbackend(:zygote)  # only AD that works right now
        using ForwardDiff
        Turing.setadbackend(:forwarddiff)   
    =#

  Logging.disable_logging(Logging.Debug-2000)  # force re-enable logging
  
  #  Run sampler, collect results.
  n_sample_test = 30
  n_adapts_test = 30
  n_chains_test = 4

  #=
      # define sampler 
      # test with SGLD as it is fast and provides behavioural range
      # ensure basic solutions are within range .. testing balance of parameter effects (positive, follows data, etc)
      # then choose NUTS .. it is probably the simplest choice:
      turing_sampler_test = Turing.SGLD()   # Stochastic Gradient Langevin Dynamics (SGLD)
      turing_sampler_test = Turing.HMC(0.01, 7)
      turing_sampler_test = Turing.SMC()
      turing_sampler_test = Turing.HMCDA(0.25, 0.65)  #  total leapfrog length, target accept ratio.
      turing_sampler_test = Turing.NUTS{Turing.ReverseDiffAD{true}}( n_adapts_test, 0.65 ) # , init_ϵ=0.001
      turing_sampler_test = Turing.NUTS( 0.65 ) # , init_ϵ=0.001
  =#

  turing_sampler_test = Turing.NUTS(n_adapts_test, 0.65; max_depth=7, init_ϵ=0.0125 )
 
  #=
      # collect good seeds. criteria: (g10

      # must re-run if redefining params:
      fmod = size_structured_dde_turing( PM=PM, solver_params=solver_params )
      
      fmod = logistic_discrete_turing_historical( PM )  # q only
 
    =#

  
  res = sample( fmod, turing_sampler_test, n_sample_test  ) # to see progress -- about 5 min
  
  #=
      res = sample( fmod, turing_sampler_test, MCMCThreads(), n_sample_test, n_chains ) 
      res = sample( fmod, turing_sampler_test, MCMCThreads(), n_sample_test, init_params=summarize(res).nt[2]  ) # test to see if restart works well
      res = sample( fmod, turing_sampler_test, n_sample_test, init_params=summarize(res).nt[2]  ) # test to see if restart works well
  =#

  showall( summarize( res ) )

  # n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
  m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(res  )

  # fishing (kt), relative Fishing mortlaity, instantaneous fishing mortality:
  Fkt, FR, FM = fishery_model_mortality() 
  pl = fishery_model_plot( toplot=("survey", "fishing", "nofishing", "trace") )
 
  #=
      describe(res)
      plot(res)
      summarystats(res)

      pl = density!(res[ "K[1]" ])  
      pl = plots_diagnostic( res, toplot="K[1]" )  # same thing 
      pl = plots_diagnostic( res, vn="K", i=1 )  # same thing 
      # savefig(pl, joinpath( model_outdir, string("diagnostic", aulab, vn, ".pdf") )  )

      summarystats(res[:,1:4,:])
      plot(res)
      plot( traceplot(res) )
      plot( meanplot(res) )
      plot( density(res) )
      plot( histogram(res) )
      plot( mixeddensity(res) )
      plot( autocorplot(res) )
      # labels = [:b[1], :b[2]]; corner(res, :b)
      #  pl = plot(pl, ylim=(0, 0.65))

      # diagnostic plots
      pl = fishery_model_plot( toplot=("survey", "fishing", "nofishing", "trace") )
      pl = fishery_model_plot( toplot="fishing_mortality" )
      pl = fishery_model_plot( toplot="harvest_control_rule" )  # hcr with fishing mortality
      pl = fishery_model_plot( toplot="trace_footprint" )
      pl = fishery_model_plot( toplot="fishing_mortality_vs_footprint" )
      pl = fishery_model_plot( toplot="harvest_control_rule_footprint" )  # hcr with fishing footprint
      pl = fishery_model_plot( toplot="number", si=1 )  # s1 as numbers
      pl = fishery_model_plot( toplot="number", si=2 )  # s2 numbers
      pl = fishery_model_plot( toplot="number", si=6 )  # female numbers
  =#


# ------------
# Production run  

  Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
 
  #=
      # params tweaks if required: most are defined in environments 
      n_adapts=200
      n_samples=500
      n_chains=4
      rejection_rate = 0.65
      max_depth = 7
      init_ϵ = 0.0125
      turing_sampler = Turing.NUTS(n_samples, rejection_rate; max_depth=max_depth, init_ϵ=init_ϵ )
  =#
 
  # choose one:
  res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains ) 

  #= alternative: restart from pre-determined means
    res  =  sample( fmod, turing_sampler, MCMCThreads(), 
      n_samples, n_chains, init_params=FillArrays.Fill(summarize(res).nt[2], n_chains) #, thinning=thinning, discard_initial=discard_initial  
    )
    #, thinning=thinning, discard_initial=discard_initial  
    # init_params=FillArrays.Fill(summarize(res).nt[2], n_chains) 
  =#

  showall( summarize( res ) )

  # save results to (model_outdir) as a hdf5  .. can also read back in R as:  h5read( res_fn, "res" )
  # directory location is created in environment
  # res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )   #defined in *.environment.jl
  @save res_fn res

    #=  to reload a save file:

      @load res_fn res

    =#
   
  # extract values into main memory:
  # n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
  m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(res ) 

  # fishing (kt), relative Fishing mortality, instantaneous fishing mortality:
  Fkt, FR, FM = fishery_model_mortality() 


  # save a few data files as semicolon-delimited CSV's for use outside Julia
  summary_fn = joinpath( model_outdir, string("results_turing", "_", aulab, "_summary", ".csv" ) )  
  CSV.write( summary_fn,  summarize( res ), delim=";" )  # use semicolon as , also used in parm names
  
  if  occursin( r"logistic_discrete", model_variation ) 
    bio_fn1 = joinpath( model_outdir, string("results_turing", "_", aulab, "_bio_fishing", ".csv" ) )  
    CSV.write( bio_fn1,  DataFrame(bio[:,:], :auto), delim=";" )  # use semicolon as , also used in parm names
  end

  if  occursin( r"size_structured", model_variation ) 
    bio_fn1 = joinpath( model_outdir, string("results_turing", "_", aulab, "_bio_fishing", ".csv" ) )  
    CSV.write( bio_fn1,  DataFrame(bio[:,:,1], :auto), delim=";" )  # use semicolon as , also used in parm names

    bio_fn2 = joinpath( model_outdir, string("results_turing", "_", aulab, "_bio_nofishing", ".csv" ) )  
    CSV.write( bio_fn2,  DataFrame( bio[:,:,2], :auto), delim=";" )  # use semicolon as , also used in parm names
  end

  fm_fn = joinpath( model_outdir, string("results_turing", "_", aulab, "_fm", ".csv" ) )  
  CSV.write( fm_fn,  DataFrame( FM, :auto), delim=";" )  # use semicolon as , also used in parm names


  # plots 
  #  pl = plot(pl, ylim=(0, 0.65))

  if  occursin( r"logistic_discrete", model_variation ) 
  
    # annual snapshots of biomass (kt) 
    pl = fishery_model_plot( toplot=("survey", "fishing" ) )
    savefig(pl, joinpath( model_outdir, string("plot_predictions_", aulab, ".pdf") )  )

    # plot fishing mortality
    pl = fishery_model_plot( toplot="fishing_mortality" )
    pl = plot(pl, ylim=(0, 0.65))
    savefig(pl, joinpath( model_outdir, string("plot_fishing_mortality_", aulab, ".pdf") )  )

    # HCR plot
    pl = fishery_model_plot( toplot="harvest_control_rule", n_sample=1000 ) #, alphav=0.01 )  # hcr
    pl = plot(pl, ylim=(0, 0.65))
    savefig(pl, joinpath( model_outdir, string("plot_hcr_", aulab, ".pdf") )  )

  end


  if  occursin( r"size_structured", model_variation ) 
    
    pl = fishery_model_plot( toplot=("survey", "fishing", "nofishing", "trace"), alphav=0.025 )
    savefig(pl, joinpath( model_outdir, string("plot_predictions_everything_", aulab, ".pdf") )  )


    pl = fishery_model_plot( toplot=("trace", "survey"), alphav=0.025 )
    savefig(pl, joinpath( model_outdir, string("plot_predictions_trace_", aulab, ".pdf") )  )

    # annual snapshots of biomass (kt) 
    pl = fishery_model_plot( toplot=("survey", "fishing" ) )
    savefig(pl, joinpath( model_outdir, string("plot_predictions_", aulab, ".pdf") )  )

    # annual snapshots of biomass (kt) 
    pl = fishery_model_plot( toplot=("survey", "fishing", "nofishing") )
    savefig(pl, joinpath( model_outdir, string("plot_predictions_full_", aulab, ".pdf") )  )

    # plot fishing mortality
    pl = fishery_model_plot( toplot="fishing_mortality" )
    pl = plot!(pl; ylim=(0, 1.4 ) )
    savefig(pl, joinpath( model_outdir, string("plot_fishing_mortality_", aulab, ".pdf") )  )

    # HCR plot
    pl = fishery_model_plot( toplot="harvest_control_rule", n_sample=1000  )  # hcr
    pl = plot!(pl; ylim=(0, 0.65 ) )
    savefig(pl, joinpath( model_outdir, string("plot_hcr_", aulab, ".pdf") )  )

    # HCR footprint
    pl = fishery_model_plot( toplot="harvest_control_rule_footprint", n_sample=1000  )  # hcr with fishing footprint
    savefig(pl, joinpath( model_outdir, string("plot_hcr_footprint_", aulab, ".pdf") )  )
        
    # fishery footprint
    pl = fishery_model_plot( toplot="footprint" )
    savefig(pl, joinpath( model_outdir, string("plot_footprint_", aulab, ".pdf") )  )

    pl = fishery_model_plot( toplot="trace_footprint", alphav=0.02 )
    savefig(pl, joinpath( model_outdir, string("plot_footprint_trace_", aulab, ".pdf") )  )
    
    pl = fishery_model_plot( toplot="fishing_mortality_vs_footprint" )
    savefig(pl, joinpath( model_outdir, string("plot_fishing_mortality_vs_footprint_", aulab, ".pdf") )  )
 
    # timeseries of predictions (number; kn and pl =plot) -- not relevent if only 1 state varable
    statevar = 1  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )

    statevar = 2  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )

    statevar = 3  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )

    statevar = 4  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )

    statevar = 5  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )
    
    statevar = 6  # index of S
    pl = fishery_model_plot( toplot="number", si=statevar )  # s1 as numbers
    savefig(pl, joinpath( model_outdir, string("plot_predictions_timeseries_",  aulab, "_", statevar, ".pdf") )  )
 
  end
 
  
  # projections relative to status quo .. update the values below 
  status_quo_tac = [0.9799, 7.345, 0.125 ][ki] # 2022
 

  for frac in (0.8, 0.9, 1.0, 1.1, 1.2)

    Catch = status_quo_tac * frac 
    m, num, bio, trace, trace_bio, trace_time = project_with_constant_catch( res, solver_params=solver_params, PM=PM, Catch=Catch, ny_fishing_pattern=5  )
    Fkt, FR, FM = fishery_model_mortality() 
    pl = fishery_model_plot( toplot=("trace_projections"), alphav=0.025, time_range_predictions=(2002,year_assessment+4)  )

    savefig(pl, joinpath( model_outdir, string("plot_trace_projections_", aulab, "__", frac, "__",".pdf") )  )
    pl = fishery_model_plot( toplot=("trace_footprint_projections"), alphav=0.025,  time_range_predictions=(2002,year_assessment+4)  )
    savefig(pl, joinpath( model_outdir, string("plot_trace_footprint_projections_", aulab, "__", frac, "__",".pdf") )  )

  end


### end
### -------------






### -------------
## FOLLOWING are tests only  of other methods of parameter estimation  (ignore)

restart_method = false
if restart_method
  using Optim, StatsBase
  res_map = Optim.optimize( fmod, MAP())  # find starting point (modes from maximum aposteriori)
  coeftable( res_map)
  res = sample( fmod, Turing.NUTS(), 30, init_params=res_map.values.array )
end



do_variational_inference = false
if do_variational_inference
  # to do Variational Inference (an sd term goes less than 0 .. not sure of the cause ):

  res_vi =  vi(fmod, Turing.ADVI( 10, 1000));

     # Run turing_sampler, collect results. @doc(Variational.ADVI) :
     # samples_per_step::Int64
     # Number of samples used to estimate the ELBO in each optimization step.
     # max_iters::Int64
     # Maximum number of gradient steps.

  res_vi_samples = rand( res_vi, 1000)  # sample via simulation


  p1 = histogram(res_vi_samples[1, :]; bins=100, normed=true, alpha=0.2, color=:blue, label="")
  density!(res_vi_samples[1, :]; label="s (ADVI)", color=:blue, linewidth=2)
  density!(res, :s; label="s (NUTS)", color=:green, linewidth=2)
  vline!([var(x)]; label="s (data)", color=:black)
  vline!([mean(res_vi_samples[1, :])]; color=:blue, label="")

  p2 = histogram(res_vi_samples[2, :]; bins=100, normed=true, alpha=0.2, color=:blue, label="")
  density!(res_vi_samples[2, :]; label="m (ADVI)", color=:blue, linewidth=2)
  density!(res, :m; label="m (NUTS)", color=:green, linewidth=2)
  vline!([mean(x)]; color=:black, label="m (data)")
  vline!([mean(res_vi_samples[2, :])]; color=:blue, label="")

  plot(p1, p2; layout=(2, 1), size=(900, 500))


  do_maximum_likelihood = false
  if do_maximum_likelihood
    res_mle = Turing.optimize(fmod, MLE())
    res_mle = Turing.optimize(fmod, MLE())
    res_mle = Turing.optimize(fmod, MLE(), NelderMead())
    res_mle = Turing.optimize(fmod, MLE(), SimulatedAnnealing())
    res_mle = Turing.optimize(fmod, MLE(), ParticleSwarm())
    res_mle = Turing.optimize(fmod, MLE(), Newton())
    res_mle = Turing.optimize(fmod, MLE(), AcceleratedGradientDescent())
    res_mle = Turing.optimize(fmod, MLE(), Newton(), Optim.Options(iterations=10_000, allow_f_increases=true))

    using StatsBase
    coeftable(res_mle)
end


do_maximum_aposteriori = false
if do_maximum_aposteriori
  res_map = Turing.optimize(fmod, MAP())
  using StatsBase
  coeftable(res_map)
end


 
rng = MersenneTwister(26)
resp = predict(rng, textmodel_marginal_pred(data), res)
@df resp ecdfplot(:"b[1]"; label="birth rate 1")



using Flux, DiffEqFlux
params = Flux.params(p)

msol =  solve( prob,  solver, callback=cb, saveat=dt)    # isoutofdomain=(y,p,t)->any(x->x<0,y)  # to force positive

function predict_rd() # Our 1-layer "neural network"
  solve(prob, solver, p=p, saveat=dt)[1] # override with new parameters  # isoutofdomain=(y,p,t)->any(x->x<0,y)  # to force positive
end

loss_rd() = sum(abs2,x-1 for x in predict_rd()) # loss function (squared absolute) x-1

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cbflux = function () #callback
  # function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=p), solver,saveat=dt), ylim=(0,kmu*2)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, params, data, opt, cb = cbflux)

m = Chain(
  Conv((2,2), 1=>16, relu),
  x -> maxpool(x, (2,2)),
  Conv((2,2), 16=>8, relu),
  x -> maxpool(x, (2,2)),
  x -> reshape(x, :, size(x, 4)),
  x -> solve(prob, solver,u0=x,saveat=0.1)[1,:],
  Dense(288, 10), softmax) |> gpu





      
 