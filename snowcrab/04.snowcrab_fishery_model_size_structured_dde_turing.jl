
# ------------------------------
# DDE size structured (Delay Differentail EQs size structured)
# https://diffeq.sciml.ai/stable/tutorials/dde_example/


if false
  # ==== R-code ====
    # prep data
    # NOTE::: this requires 03.snowcrab_carstm.r to be completed
    source( file.path( code_root, "bio_startup.R" )  )
    loadfunctions("bio.snowcrab")
    # fishery landings has a weekly time step = 2/52 ~ 0.0385 ~ 0.04  X dt=0.01 seems to work best
    fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics", for_julia=TRUE, time_resolution=2/52  )
  # ==== R-code ====
end


if false
  # in REPL< this needs to be loaded first as it skips the startup.jl
  project_directory = @__DIR__() #  same folder as the current file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
end


# ---------------
# run-level options

  # choose a region of interest"
  aulab ="cfanorth"
  aulab ="cfasouth"
  aulab ="cfa4x"


  yrs = 1999:2021  # <<<<<<<<-- change

# load libs and options and prepare data for diffeq/turing model and set default parameters
  include( "size_structured_dde_environment.jl" )

# to reload dynamical core model, turing estimation model, etc
#  include( "size_structured_functions.jl" )

# to test dynamical model with generic/random parameters:
# include(  "size_structured_dde_test.jl" )


# run model settings / options / overrides to prep for DifferentialEquations and Turing

  # solver = MethodOfSteps(Tsit5())  # usually works well. Alt:  solver = MethodOfSteps(Rodas5())

  if false
    ## test DifferentialEquations DDE model -- ignore
    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    b=[2.0, 1.0]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;

    d=[0.1, 0.2, 0.2, 0.2, 0.2, 0.2];
    d2=[0.1, 0.2, 0.2, 0.2, 0.2, 0.2];
    v=[0.8, 0.9, 0.9, 0.9];
    u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ]  ;
   # u0 = [0.796203667247763, 0.5370253706021624, 0.49999994680822296, 0.4999999261197093, 0.3518981144352651, 0.7221526754207315]
    params = ( b, K, d, d2, v,  tau, hsa)
    prob = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]
    msol2 =  solve( prob,  solver, callback=cb, saveat=dt )
    plot( msol2, ; legend=false, xlim=(1999,2021), label="test" )
  end



  if false
    ## test Turing model and timings  -- ignore
    Random.seed!(1)

    n_samples = 20
    n_adapts = 20

    res  =  sample( fmod, MH(), 100  )  #  Metropolis-Hastings
    res  =  sample( fmod, DynamicNUTS(), n_samples  )

    leapfrog_stepsize = 0.01
    n_leapfrog_steps = 50
    res  =  sample( fmod, Turing.HMC(leapfrog_stepsize, n_leapfrog_steps ), n_samples  )

    res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65 ), n_samples  )

    # include( "size_structured_dde_environment.jl" )
    res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65, init_ϵ=0.05, max_depth=7), n_samples,  progress=true, drop_warmup=true  )

    showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

  end

# production  uses Turing.NUTS
Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info

n_samples = 1000  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs;; 29hrs 100/100 cfasouth
n_adapts = 1000
n_chains = 5  #   # n_chains = Threads.nthreads()
turing_sampler = Turing.NUTS(n_adapts, 0.65; max_depth=7, init_ϵ=0.05)  ;# stepsize based upon previous experience
res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains )
# if on windows and threads are not working, use single processor mode:
# res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

# save results as a hdf5
 fn = joinpath( project_directory, "ignore", string("size_structured_dde_turing_data", "_", aulab, ".hdf5" ) )
 @save fn res
# @load fn res

# can read back in R as:
# h5read( paste("/home/jae/julia/snowcrab/ignore/size_structured_dde_turing_data", "_", aulab, ".hdf5"), "res")


# ------------------------------
# plots

  plot()
  vn = "model_sd"; density!(res[ Symbol(vn) ])


  plot()
  vn = "b[1]"; density!(res[ Symbol(vn) ])
  vn = "b[2]"; density!(res[ Symbol(vn) ])

  plot()
  vn = "d[1]"; density!(res[ Symbol(vn) ])
  vn = "d[2]"; density!(res[ Symbol(vn) ])
  vn = "d[3]"; density!(res[ Symbol(vn) ])
  vn = "d[4]"; density!(res[ Symbol(vn) ])
  vn = "d[5]"; density!(res[ Symbol(vn) ])
  vn = "d[6]"; density!(res[ Symbol(vn) ])

  plot()
  vn = "K[1]"; density!(res[ Symbol(vn) ])
  vn = "K[2]"; density!(res[ Symbol(vn) ])
  vn = "K[3]"; density!(res[ Symbol(vn) ])
  vn = "K[4]"; density!(res[ Symbol(vn) ])
  vn = "K[5]"; density!(res[ Symbol(vn) ])
  vn = "K[6]"; density!(res[ Symbol(vn) ])


# plot simulation traces of FB with and without fishing
  (trace_nofishing, trace_fishing, pl) = size_structured_predictions( res; ns=500, plot_k=1, alpha=0.1 )  # model traces
  (pl)

  savefig(pl, joinpath( project_directory, "ignore", string("size_structured_dde_turing_plots_sim_", aulab, ".pdf") )  )


  # annual snapshots of numerical abundance (relative number) converted to biomass (kt)
  (m, num, bio)  = size_structured_predictions_annual(res; prediction_time=prediction_time, ns=500)

  # extract sims (with fishing)
  g = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
  size(g)

  # plot biomass
  # gr()
  pl = plot()
  pl = plot!(pl, prediction_time, g;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  k = 1
  yhat = ( S[:,k] .- mean(res[:,Symbol("qc[$k]"),:]  ) ) ./ mean(res[:,Symbol("q[$k]"),:]) .* mean(res[:,Symbol("K[$k]"),:]  )

  # S[i,k] ~ TruncatedNormal( msol.u[ii][k] * q[k] + qc[k], bpsd, 0.0, 1.0)

  if nameof(typeof(mw)) == :ScaledInterpolation
    yhat = yhat .* mw(yrs) ./ 1000.0  ./ 1000.0
  else
    yhat = yhat .* scale_factor
  end
  pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )

  savefig(pl, joinpath( project_directory, "ignore", string("size_structured_dde_turing_plots_predictions_", aulab, ".pdf") )  )



  # plot numbers
  # gr()
  k = 5
  gk = num[:,k,:,1]
  pl = plot()
  pl = plot!(pl, prediction_time, gk;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  yhat = ( S[:,k] .- mean(res[:,Symbol("qc[$k]"),:]  ) ) ./ mean(res[:,Symbol("q[$k]"),:]) .* mean(res[:,Symbol("K[$k]"),:]  )

  pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )




  # plot fishing mortality
  removed_annual = removals_aggregate( removed, fish_year )

  Fkt = removed_annual[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  FR =  Fkt ./ ( Fkt .+  g[1:length(survey_time),:] )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F

  pl = plot()
  pl = plot!(pl, survey_time, FM ;  alpha=0.1, color=:lightslateblue)
  pl = plot!(pl, survey_time, mean(FM, dims=2) ;  alpha=0.8, color=:slateblue, lw=4)
  pl = plot!(pl, xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  pl = plot!(pl, ylim=(0, maximum(FM)*1.1 ) )
  pl = plot!(pl ; legend=false )

  savefig(pl, joinpath( project_directory, "ignore", string("size_structured_dde_turing_plots_fishingmortality_", aulab, ".pdf") )  )



  # HCR plot
  pl = plot()

  # mean weight by year
  sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor

  # sample and plot posterior K
  K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass

  nsample = 500
  o = rand(K, nsample)
  pl = vline!(pl, o;  alpha=0.05, color=:limegreen )
  pl = vline!(pl, o./2;  alpha=0.05, color=:darkkhaki )
  pl = vline!(pl, o./4;  alpha=0.05, color=:darkred )

  pl = vline!(pl, [mean(o)];  alpha=0.6, color=:chartreuse4, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

  pl = vline!(pl, [mean(o)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )

  pl = vline!(pl, [mean(o)/4.0];  alpha=0.6, color=:darkred, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )


  nt = length(survey_time)
  ns = size(g)[2]
  gs = g[1:nt,:]

  colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]

  # scatter!( gs, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)

  gg = mean(gs, dims=2)
  ff = mean(FM, dims=2)

  # scatter!( [gs[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
  pl = plot!(pl, gg, ff ;  alpha=0.8, color=:slateblue, lw=3)

  pl = scatter!(pl,  gg, ff ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
    series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
  pl = scatter!(pl,  [gg[nt]], [ff[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
  ub = max( quantile(o, 0.95), maximum( gg ) ) * 1.05
  pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(ff ) * 1.05  ) )

  # add predictions ???
  gs = g[nt,:]

  sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(prediction_time) ./ 1000.0  ./ 1000.0 : scale_factor

  savefig(pl, joinpath( project_directory, "ignore", string("size_structured_dde_turing_plots_hcr_", aulab, ".pdf") )  )




  # parameter estimates for output
  MSY    = r * exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY




# rng = MersenneTwister(26)
# resp = predict(rng, textmodel_marginal_pred(data), res)
# @df resp ecdfplot(:"b[1]"; label="birth rate 1")



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
end







using Flux, DiffEqFlux
params = Flux.params(p)

msol =  solve( prob,  solver, callback=cb, saveat=dt)    # isoutofdomain=(y,p,t)->any(x->x<0,y)  # to force positive

function predict_rd() # Our 1-layer "neural network"
  solve(prob, solver,p=p,saveat=dt)[1] # override with new parameters  # isoutofdomain=(y,p,t)->any(x->x<0,y)  # to force positive
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







if false
  # should the models be stored in files, running this way is one option
  # discrete version
  include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_basic.jl")

  # ode version
  include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_ode.jl")
end



# more inits
Random.seed!(1);

test = fishery_model_turing( Y, Kmu, Ksd, removals)

# sample from it as a test and precompile rquired functions
res_smc =  sample( fishery_model_turing( Y, Kmu, Ksd, removals ),  SMC(), 10 )

do_variational_inference = false
if do_variational_inference
  # to do Variational Inference (an sd term goes less than 0 .. not sure of the cause ):

  using Flux, Turing
  res_vi =  vi(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.ADVI(10, 1000));

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

end


do_maximum_likelihood = false
if do_maximum_likelihood
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), NelderMead())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), SimulatedAnnealing())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), ParticleSwarm())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), Newton())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), AcceleratedGradientDescent())
    res_mle = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MLE(), Newton(), Optim.Options(iterations=10_000, allow_f_increases=true))

    using StatsBase
    coeftable(res_mle)
end


do_maximum_aposteriori = false
if do_maximum_aposteriori
  res_map = optimize(fishery_model_turing( Y, Kmu, Ksd, removals ), MAP())
end


do_mcmc = false
if do_mcmc

    # Turing.setadbackend(:reversediff)
    # Turing.setadbackend(:forwarddiff)  # slightly faster .. (default)
    # Turing.setadbackend(:zygote)  #  does not like loops
    # Turing.setadbackend(:tracker) # does not like loops

    # Load Distributed to add processes and the @everywhere macro.
    # addprocs(4)

     # Initialize everything on all the processes.
     # Note: Make sure to do this after youve already loaded Turing,
     #     so each process does not have to precompile.
     #     Parallel sampling may fail silently if you do not do this.
     @everywhere using Turing


    # check no threads
    Threads.nthreads()

    # SMC: number of particles.
    # PG: number of particles, number of iterations.
    # HMC: leapfrog step size, leapfrog step numbers.
    # Gibbs: component turing_sampler 1, component turing_sampler 2, ...
    # HMCDA: total leapfrog length, target accept ratio.
    # NUTS: number of adaptation steps (optional), target accept ratio.

    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.SMC(), 5000 )
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.SMC(), MCMCThreads(), 1000, 4)
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.PG(10), 1000)
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.HMC(0.1, 5), 1000)
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.Gibbs(PG(10, :m), HMC(0.1, 5, :s²)), 1000)
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.HMCDA(0.15, 0.65), 1000)
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.NUTS(0.65), 1000)

    # multithreaded:
    res = sample(fishery_model_turing( Y, Kmu, Ksd, removals ), Turing.NUTS(0.65), MCMCThreads(), 1000, 4 )  # 30 min

    warmup = 1000
    niters = 4000

    res =  sample( fishery_model_turing( Y, Kmu, Ksd, removals ),  Turing.SMC(),   niters, n_adapts=warmup, progress=false, verbose=false, drop_warmup=true )

    # res =  sample( fishery_model_turing( Y, Kmu, Ksd, removals ),  PG(10L), niters, n_adapts=warmup, progress=false, verbose=false, drop_warmup=true )
    res =  sample( fishery_model_turing( Y, Kmu, Ksd, removals ),  Turing.NUTS(0.65),   niters, n_adapts=warmup, progress=false, verbose=false, drop_warmup=true )

    describe(res)

    p = plot(res)


    summarize( res)

    summarystats(res[:,1:4,:])

    using  Plots, LaTeXStrings, StatsPlots

    plot(res)


    plot( traceplot(res) )
    plot( meanplot(res) )
    plot( density(res) )
    plot( histogram(res) )
    plot( mixeddensity(res) )
    plot( autocorplot(res) )


    labels = [:r, :K, :q, :qc]
    corner(res, labels)


end



# Calculate MSE for our test set.
loss = sum((predictions - test_label) .^ 2) / length(test_label)


