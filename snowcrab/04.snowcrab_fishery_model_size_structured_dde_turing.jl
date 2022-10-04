 
# ------------------------------
# DDE size structured (Delay Differentail EQs size structured)
# https://diffeq.sciml.ai/stable/tutorials/dde_example/


# ------------------------------
# run-level options

  nS = 6  # no. state variables

  yrs = 1999:2021  # <<<<<<<<-- change

  nT = length(yrs)
  nP = 5  # number of predictions into future (with no fishing)
  nM = nP + nT  # total number of prediction years

  dt = 0.01 # time resolution of diff eq model solutions

    
  # choose a region of interest"

  aulab ="cfanorth"

  aulab ="cfasouth"

  aulab ="cfa4x"
      

{
  # ==== R-code to prep data ======
  
  # NOTE::: require 03.snowcrab_carstm.r to be completed 
 
  source( file.path( code_root, "bio_startup.R" )  )
  loadfunctions("bio.snowcrab")
  fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/52)

}
        

if false
  # if doing manual startup
  project_directory = @__DIR__() #  same folder as the file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
end

# ------------------------------
# load libs and check settings

# add Turing@v0.21.10
 
pkgs = [ 
  "Revise", "MKL", "Logging", "StatsBase", "Statistics", "Distributions", "Random",   
  # "DynamicHMC", "AdvancedHMC",  "AdvancedMH",  "DynamicPPL",  "AbstractPPL",  "Memoization", 
  "ForwardDiff",
  "Plots", "StatsPlots", "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
  "Turing", "ModelingToolkit", "DifferentialEquations", "Interpolations", "LinearAlgebra" 
]

for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages
     

 
# prepare data for diffeq/turing model and set default parameters
include( "size_structured_dde_turing_data.jl" )
 
# load dynamical core model
include( "size_structured_dde!.jl" ) 

# to test dynamical model with generic parameters:
# include( joinpath( project_directory, "size_structured_dde_test.jl" ))    

# turing estimation model 
include( "size_structured_dde_turing.jl" ) 


# ---------------
# run model settings / options / overrides
  
solver = MethodOfSteps(Tsit5())  
  # solver = MethodOfSteps(Rodas5())   # default

prob = DDEProblem( size_structured_dde!, u0, h, tspan, p, constant_lags=lags  )
fmod = size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM, solver, dt )

 
if false
    # for testing and timings
    n_samples = 12
    n_adapts = 12

    Random.seed!(1)
    
    res  =  sample( fmod, MH(), 100  )  #  Metropolis-Hastings 
    
    res  =  sample( fmod, DynamicNUTS(), n_samples  )
 
    leapfrog_stepsize = 0.01
    n_leapfrog_steps = 50 
    res  =  sample( fmod, Turing.HMC(leapfrog_stepsize, n_leapfrog_steps ), n_samples  )
 
    res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65 ), n_samples  )
    
    res  =  sample( fmod, Turing.NUTS(n_adapts, 0.65, init_ϵ=0.025, max_depth=7), n_samples,  progress=true, drop_warmup=true  )
     
    show(stdout, "text/plain", summarize(res)) # display all estimates

end


# production  

Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info

n_samples = 1000  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs;; 29hrs 100/100 cfasouth
n_adapts = 1000
n_chains = 4
# n_chains = Threads.nthreads() 
 
turing_sampler = Turing.NUTS(n_adapts, 0.65; max_depth=7, init_ϵ=0.05)  ;# stepsize based upon previous experience

res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains )
# if on windows and threads are not working, use single processor mode:
# res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)


# ------------------------------
# save results as a hdf5
using JLD2

fn = joinpath( project_directory, string("size_structured_dde_turing_data", "_", aulab, ".hdf5" ) )
@save fn res

# @load fn res
# can read back in R as:  
# h5read( paste("/home/jae/julia/snowcrab/size_structured_dde_turing_data", "_", aulab, ".hdf5"), "res")

# fnk = joinpath( project_directory, string("size_structured_dde_turing_data", "_", aulab, "_K_", ".hdf5" ) )
# tt =res[:,Symbol("K[1]"),:]
# @save fnk tt

# ------------------------------
# process outputs
# display all estimates
show(stdout, "text/plain", summarize(res))




# -------------------------
# plot timeseries of mean fields
include( "size_structured_dde_turing_plot.jl" )
 

plot() 
 
# pdensity = size_structured_dde_turing_plot( selection="posterior_density", vn="b[1]"  )  

plots_sim = size_structured_dde_turing_plot( selection="S K predictions predictionmeans", si=[1], mw=mw ) 
plot!(; xlim=(minimum(yrs)-1.5, maximum(yrs)+7.5  ) )

# display(plots_sim)
savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".png") )
# savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".pdf") )
# savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".svg") )


plot(0)
plots_fishing = size_structured_dde_turing_plot( selection="S K withfishing withoutfishing ", si=[1], mw=mw)
# display(plots_fishing)
savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".png") ) 
# savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".pdf") ) 
# savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".svg") ) 



# ------------------------------
# misc computed quantities

chains_params = Turing.MCMCChains.get_sections(res, :parameters)
generated_quantities(fmod2, chains_params)

o = generated_quantities(fmod, res) 
  
pm = ( b=b, K=K, d=d, v=v, tau=tau, hsa=hsa ) 


function fishing_mortality( removed, abundance )
    -log(  1.0 - (removed  / abundance)  )  ;
end

function size_structured_predictions( res)
 
  nchains = size(res)[3]
  nsims = size(res)[1]
  
  m = zeros(nM, nS, 2, nchains*nsims)
 
  z = 0

  gr() 
  theme(:default)

  pl =plot()
  
  for j in 1:nsims  # nsims 
  for l in 1:nchains #nchains
    z += 1
  #  for k in 1:nS  # nstate vars

        b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
        K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
        d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
        v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]

        q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
        qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]

        u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

        pm = ( b, K, d, v, tau, hsa ) 
        
        prb = remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm ) 
        
        msol = solve( prb, solver, callback=cb, saveat=dt ) 
        msol2 = solve( prb, solver, saveat=dt ) # no call backs

        yval = reduce(hcat, msol.u)'[:,1]
        yval2 = reduce(hcat, msol2.u)'[:,1]

        if nameof(typeof(mw)) == :ScaledInterpolation
            yval = yval .* mw(msol.t) ./ 1000.0  ./ 1000.0 
            yval2 = yval2 .* mw(msol2.t) ./ 1000.0  ./ 1000.0 
        else
            yval = yval .* scale_factor
            yval2 = yval2 .* scale_factor
        end

        yval = vec(yval)
        yval2 = vec(yval2)

        pl = plot!( pl, msol.t, yval, alpha=0.1, lw=2, color=:teal ) 
        pl = plot!( pl, msol2.t, yval2, alpha=0.1, lw=2, color=:lightgreen ) 
     
    
        for i in 1:nM
            ii = findall(x->x==prediction_time[i], msol.t)[1]
            jj = findall(x->x==prediction_time[i], msol2.t)[1]
            m[i,:,1,z] = msol.u[ii]    
            m[i,:,2,z] = msol2.u[jj]
        end

    #  end
             

  end
  end
    pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
    pl =  plot!(pl; ylim=(0, maximum(yval)*1.1 ) )
    pl =  plot!(pl; legend=false )

  return(m,pl)
  
 end    
   
m, pl = size_structured_predictions( res)






g[:,1,z] = m[:,k,z]
g[:,2,z] = m2[:,k,z]

if nameof(typeof(mw)) == :ScaledInterpolation
    g[:,1,z] .*=  mw(prediction_time) ./ 1000.0  ./ 1000.0 
    g[:,2,z] .*=  mw(prediction_time) ./ 1000.0  ./ 1000.0 
else
    g[:,1,z] .*=  scale_factor
    g[:,2,z] .*=  scale_factor
end

# plot!(prediction_time, g[:,1,z];  alpha=0.2, color=:orange)
# plot!(prediction_time, g[:,2,z];  alpha=0.2, color=:orange)





for i in 1:length(res)  
    w = zeros(nM)
    for j in 1:nM
        w[j] = res[i, Symbol("K[$k]"),l] * res[i, Symbol("m[$j,$k]"),l] 
    end
    if nameof(typeof(mw)) == :ScaledInterpolation
        w = w .* mw(prediction_time) ./ 1000.0  ./ 1000.0 
    else
        w = w .* scale_factor
    end
    plot!(prediction_time, w;  alpha=0.1, color=:orange)
end


 
  # params need to be named  .. return only that which is specified by "return()", below
  # deterministic computations: do from similations:
 
  FM = zeros(nT+M)
  BM = zeros(nT+M)
  # CA = zeros(nT+M)

  # CA[1:nT] = removed ./ K
  # CA[(nT+1):(M+nT)] = er .* bm[(nT):(M+nT-1)]
  # CA = 1.0 .- CA / bm

  F =  -log( max.(C, smallnumber) )  ;
  
  # parameter estimates for output
  MSY    = r * exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
   
 
    m = Matrix{T}(undef, nM, nS)   # 
   
    pm = ( b, K, d, v, tau, hsa ) 
    # @show a

    # process model
    msol = solve( 
        remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm ), 
        solver, 
        callback=cb, 
        # isoutofdomain=(y,p,t)->any(x->x<0, y), 
        saveat=dt
    ) 
    
    # @show msol
    
    if msol.retcode != :Success 
      Turing.@addlogprob! -Inf
      return nothing
    end
   
    # for i in 1:nM
    #     ii = findfirst(x->isapprox(prediction_time[i], x), msol.t)[1]
    #     for k in 1:nS
    #         m[i,k] ~ TruncatedNormal( msol.u[ii][k], bpsd, 0.01, 0.99 )  
    #     end
    # end
  end
   
gq = generated_quantities(fmod2, res)

# gq = generated_quantities(res, values(p), keys(p))
# m = vec(getindex.( gq, 1))
# density!(m, lab="generated quantity (VI)")
# vline!([0], lab="true value")

 

# look at predictions:
M0_pred = Vector{Union{Missing, Float64}}(undef, length(S1))
fmod_pred = fmod( M0_pred, kmu,  tspan, prob  ) 
 
predictions = predict(fmod_pred, res)
y_pred = vec(mean(Array(group(predictions, :S1)); dims = 1));

plot( S1, y_pred )
sum(abs2, S1 - y_pred) ≤ 0.1


 
 
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
 
   
