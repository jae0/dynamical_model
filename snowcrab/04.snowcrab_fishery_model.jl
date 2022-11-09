
# ------------------------------
# fishery model

# run-level options

# choose model and area
# model_variation = "logistic_discrete_basic"
# model_variation = "logistic_discrete_map"
# model_variation = "logistic_discrete"  # default (for discrete)
# model_variation = "size_structured_dde_unnormalized"  # basic model without normaliztion
# model_variation = "size_structured_dde_normalized"  # default (for continuous)

model_variation = "size_structured_dde_normalized"  # default (for continuous)


# choose a region of interest"
aulab ="cfanorth"
aulab ="cfasouth"
aulab ="cfa4x"


yrs = 1999:2021  # <<<<<<<<-- change


# create data (in R)
if false

  #  https://mybinder.org/v2/gh/jae0/dynamical_model/main

    # ==== R-code ==== start
        # # NOTE::: this requires 03.snowcrab_carstm.r to be completed
        # source( file.path( code_root, "bio_startup.R" )  )
        # loadfunctions("bio.snowcrab")
        # # prep data for discrete version
        # if (grepl("logistic_discrete", model_variation )) {
        #     fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics", for_julia=TRUE   )
        # }
        # # prep data for continuous version: 
        # if (grepl("size_structured", model_variation)) {
        #     # fishery landings has a weekly time step = 2/52 ~ 0.0385 ~ 0.04  X dt=0.01 seems to work best
        #     fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics", for_julia=TRUE, time_resolution=2/52  )
        # }
    # ==== R-code ==== end
  
  # in REPL or VSCODE, this needs to be loaded first as "startup.jl" is skipped
  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )

      
      
end


# ---------------
# load libs and options and prepare data for diffeq/turing model and set default parameters
include( joinpath( project_directory, "fishery_model_environment.jl"  ))  # bootstrap different project environments depending on above choices




debugging = false
if debugging
    # if debugging/development:
     
    # to test dynamical model with generic/random parameters
    if @isdefined fishery_model_test  
      (test, pl) = fishery_model_test( "basic" )
      pl

      (test, pl) = fishery_model_test( "random_external_forcing"  )
      pl
      
      (test, pl) = fishery_model_test( "fishing"  )
      pl
      
      (test, pl) = fishery_model_test( "nofishing" )
      pl
      showall( summarize( test ) )
    end


    res  =  sample( fmod, Turing.NUTS(30, 0.65; max_depth=7, init_ϵ=0.05), 30 ) # to see progress -- about 5 min
    # res = fishery_model_inference( fmod, n_adapts=30, n_samples=30, n_chains=1, max_depth=7, init_ϵ=0.01  )
 
    (m, num, bio, pl)  = fishery_model_predictions(res; prediction_time=prediction_time, n_sample=30 )
    fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete
    (pl)
 
    # trace plot .. only useful in continuous models, otherwise identical to predictions
    (trace_nofishing, trace_fishing, pl) = fishery_model_predictions_trace( res; n_sample=50, plot_k=1, alpha=0.1, plot_only_fishing=false )  # model traces
    (pl)

    # plot fishing mortality
    (Fkt, FR, FM, pl) = fishery_model_mortality( removed, fb, n_sample=500 ) 
    pl

    # HCR plot
    (K, bi, fm, fmsy, pl) = fishery_model_harvest_control_rule(res, yrs; FM=FM, fb=fb, n_sample=500)
    pl
   
    showall( summarize( res ) )
 
    # describe(res)
    # plot(res)
    # summarystats(res)
end



# params defined in environments .. about 3-5 hrs each
res = fishery_model_inference( fmod, rejection_rate=rejection_rate, n_adapts=n_adapts, n_samples=n_samples, 
    n_chains=n_chains, max_depth=max_depth, init_ϵ=init_ϵ )


# save results to (directory_output) as a hdf5  # directory location is created in environment
# can also read back in R as:  h5read( save_fn, "res" )
save_fn = joinpath( directory_output, string("size_structured_dde_turing_data", "_", aulab, ".hdf5" ) ) 
@save save_fn res
# @load save_fn res


    
# summaries and plots 

(directory_output)  #check output directory
# some vars
vn = "model_sd"; 
vn = "K"
vn = "K[1]" 
vn = "K[6]" 
vn = "r"
vn = "b[1]" 
vn = "b[2]" 
vn = "d[1]" 
vn = "d[6]" 

pl = density!(res[ Symbol(vn) ])  
pl = plots_diagnostic( res, vn )  # same thing 
# savefig(pl, joinpath( directory_output, string("diagnostic", aulab, vn, ".pdf") )  )

summarystats(res[:,1:4,:])

plot(res)

plot( traceplot(res) )
plot( meanplot(res) )
plot( density(res) )
plot( histogram(res) )
plot( mixeddensity(res) )
plot( autocorplot(res) )

# labels = [:b[1], :b[2]]; corner(res, :b)


# annual snapshots of biomass (kt); return m=normalized abundance, num=numbers, bio=biomass and pl=plot, where possible
(m, num, bio, pl)  = fishery_model_predictions(res; prediction_time=prediction_time, n_sample=500)
fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete
pl
savefig(pl, joinpath( directory_output, string("plot_predictions_", aulab, ".pdf") )  )

# plot fishing mortality
(Fkt, FR, FM, pl) = fishery_model_mortality( removed, fb, n_sample=500 ) 
pl
savefig(pl, joinpath( directory_output, string("plot_fishing_mortality_", aulab, ".pdf") )  )


# HCR plot
(K, bi, fm, fmsy, pl) = fishery_model_harvest_control_rule(res, yrs; FM=FM, fb=fb, n_sample=500)
pl
savefig(pl, joinpath( directory_output, string("plot_hcr_", aulab, ".pdf") )  )



if occursin.( r"size_structured", model_variation )
  # plot simulation traces of FB with and without fishing .. only for continuous models, otherwise identical to predictions (below)
  (trace_nofishing, trace_fishing, pl) = fishery_model_predictions_trace( res; n_sample=50, plot_k=1, alpha=0.1, plot_only_fishing=false )  # model traces
  pl
  savefig(pl, joinpath( directory_output, string("plot_predictions_trace_", aulab, ".pdf") )  )

  # timeseries of predictions (number; kn and pl =plot) -- not relevent if only 1 state varable
  statevar = 1  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  statevar = 2  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  statevar = 3  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  statevar = 4  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  statevar = 5  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  
  statevar = 6  # index of S
  (numS, pl)  = fishery_model_predictions_timeseries(num; prediction_time=prediction_time, plot_k=statevar )
  pl
  
  savefig(pl, joinpath( directory_output, string("plot_predictions_timeseries_", aulab, statevar, ".pdf") )  )
end


 

### end
### -------------


## FOLLOWING are tests only  of other methods of parameter estimation 


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





      
 