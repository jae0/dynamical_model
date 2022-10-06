 
# ------------------------------
# DDE size structured (Delay Differentail EQs size structured)
# https://diffeq.sciml.ai/stable/tutorials/dde_example/


if false
  # ==== R-code ====
    # prep data
    # NOTE::: this requires 03.snowcrab_carstm.r to be completed 
    source( file.path( code_root, "bio_startup.R" )  )
    loadfunctions("bio.snowcrab")
    fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/52  )  
  # ==== R-code ====
end


# ---------------
# run-level options

  yrs = 1999:2021  # <<<<<<<<-- change

  nT = length(yrs)
  nP = 5  # number of predictions into future (with no fishing)
  nM = nP + nT  # total number of prediction years

  dt = 0.01 # time resolution of diff eq model solution saves (0.02 ~ every week)
  nS = 6  # no. state variables

  # choose a region of interest"
  aulab ="cfanorth"
  aulab ="cfasouth"
  aulab ="cfa4x"
      

# load libs and options
  include( "size_structured_dde_environment.jl" )
       
# prepare data for diffeq/turing model and set default parameters
  include( "size_structured_dde_turing_data.jl" )
 
# load dynamical core model, turing estimation model, etc 
  include( "size_structured_functions.jl" ) 

# to test dynamical model with generic parameters:
# include(  "size_structured_dde_test.jl" )    
 
# run model settings / options / overrides to prep for DifferentialEquations and Turing
  solver = MethodOfSteps(Tsit5())  # usually works well. Alt:  solver = MethodOfSteps(Rodas5())  
  prob = DDEProblem( size_structured_dde!, u0, h, tspan, p, constant_lags=[tau]  )
  fmod = size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM, solver, dt )
  
    if false
      # for testing timings and model viability -- ignore
      Random.seed!(1)
  
      n_samples = 12
      n_adapts = 12
      
      res  =  sample( fmod, MH(), 100  )  #  Metropolis-Hastings 
      res  =  sample( fmod, DynamicNUTS(), n_samples  )
  
      leapfrog_stepsize = 0.01
      n_leapfrog_steps = 50 
      res  =  sample( fmod, Turing.HMC(leapfrog_stepsize, n_leapfrog_steps ), n_samples  )
      res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65 ), n_samples  )
      res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65, init_ϵ=0.05, max_depth=7), n_samples,  progress=true, drop_warmup=true  )
      
      showall(res)  # show(stdout, "text/plain", summarize(res)) # display all estimates

    end

# production  

Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info

n_samples = 1000  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs;; 29hrs 100/100 cfasouth
n_adapts = 1000
n_chains = 4  #   # n_chains = Threads.nthreads() 
 
turing_sampler = Turing.NUTS(n_adapts, 0.65; max_depth=7, init_ϵ=0.05)  ;# stepsize based upon previous experience

res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains )
# if on windows and threads are not working, use single processor mode:
# res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)


showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates


# ------------------------------
# save results as a hdf5
 fn = joinpath( project_directory, "ignore", string("size_structured_dde_turing_data", "_", aulab, ".hdf5" ) )
 @save fn res
# @load fn res

# can read back in R as:  
# h5read( paste("/home/jae/julia/snowcrab/ignore/size_structured_dde_turing_data", "_", aulab, ".hdf5"), "res")



# ------------------------------
# plots 

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


  plots_sim = size_structured_dde_turing_plot( selection="S K predictions predictionmeans", si=[1], mw=mw ) 
  plot!(; xlim=(minimum(yrs)-1.5, maximum(yrs)+7.5  ) )

  # display(plots_sim)
  savefig(plots_sim, "ignore", string("size_structured_dde_turing_plots_sim_", aulab, ".png") ) 


  plot(0)
  plots_fishing = size_structured_dde_turing_plot( selection="S K withfishing withoutfishing ", si=[1], mw=mw)
  # display(plots_fishing)
  savefig(plots_fishing, "ignore",  string("size_structured_dde_turing_plots_fishing_", aulab, ".png") )  
 
    
  o = size_structured_predictions( res; n=500, k=1 )  # model traces
    
  # annual snapshots of numerical abundance (relative number) 
  m = size_structured_predictions_annual(res; prediction_time=prediction_time, n=100)

  # extract sims (with fishing)
  g = m[:,1,:,1]   # [ yr, statevariable, sim, (with fishing=1; nofishing=2) ] 
  size(g)

  # convert number to biomass (kt)
  g .*= nameof(typeof(mw)) == :ScaledInterpolation ?  mw(prediction_time) ./ 1000.0  ./ 1000.0 : scale_factor

  # plot biomass
  gr()
  plot()
  plot!( prediction_time, g;  alpha=0.1, color=:lightslateblue)
  plot!( prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  plot!(; legend=false )
  plot!(; ylim=(0, maximum(g)*1.01 ) )
  

  # plot fishing mortality
  Fkt, FR, FM = fishing_mortality(  removed, fish_time, g[1:length(survey_time),:] )
  plot()
  plot( survey_time, FM ;  alpha=0.1, color=:lightslateblue)
  plot!( survey_time, mean(FM, dims=2) ;  alpha=0.8, color=:slateblue, lw=4)
  plot!( xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
  plot!( ; legend=false )


  # HCR plot
  plot()

  # mean weight by year
  sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor  
  
  # sample and plot posterior K
  K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass

  nsample = 500
  o = rand(K, nsample)
  vline!(o;  alpha=0.05, color=:limegreen )
  vline!(o./2;  alpha=0.05, color=:darkkhaki )
  vline!(o./4;  alpha=0.05, color=:darkred )
  
  vline!([mean(o)];  alpha=0.6, color=:chartreuse4, lw=5 )
  vline!([quantile(o, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  vline!([quantile(o, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

  vline!([mean(o)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
  vline!([quantile(o, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  vline!([quantile(o, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )

  vline!([mean(o)/4.0];  alpha=0.6, color=:darkred, lw=5 )
  vline!([quantile(o, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  vline!([quantile(o, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )

 
  nt = length(survey_time)
  ns = size(g)[2]
  gs = g[1:nt,:]

  colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]

  # scatter!( gs, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
 
  gg = mean(gs, dims=2)
  ff = mean(FM, dims=2)
  
  # scatter!( [gs[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
  plot!( gg, ff ;  alpha=0.8, color=:slateblue, lw=3)
  
  scatter!( gg, ff ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0, 
    series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
  scatter!( [gg[nt]], [ff[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
  plot!( ; legend=false )
  plot!(; xlim=(0, quantile(o, 0.98) ) )
  plot!(; ylim=(0, maximum( ff ) * 1.05 ) )

  # add predictions
  gs = g[nt,:]

  sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(prediction_time) ./ 1000.0  ./ 1000.0 : scale_factor
  
  
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



  
  ----------__--




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
 
   
