 
# ------------------------------
# DDE size structured (Delay Differentail EQs size structured)
# https://diffeq.sciml.ai/stable/tutorials/dde_example/

# NOTE::: require 03.snowcrab_carstm.r to be completed 


if false
  # if doing manual startup
  project_directory = @__DIR__() #  same folder as the file
  push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
  include( "startup.jl" )
end

# ------------------------------
# load libs and check settings
pkgs = [ 
  "Revise", "MKL", "StatsBase", "Statistics",  "Distributions", "LinearAlgebra",  "Interpolations", 
  "Plots", "StatsPlots", "MultivariateStats", "RData",
  "Turing",  "ModelingToolkit", "DifferentialEquations",  
  "StaticArrays", "LazyArrays", "FillArrays",
  "ForwardDiff", "DynamicHMC", 
  "JLD2", "HDF5"
  # "DiffResults", "Memoization", "DynamicPPL", "AbstractPPL", "AdvancedHMC", "MCMCChains", "SciMLSensitivity",
   #"Tracker" #, "ReverseDiff", "Zygote", "ForwardDiff", "Diffractor", "Memoization",
]
  
for pk in pkgs; @eval using $(Symbol(pk)); end   # Pkg.add( pkgs ) # add required packages
 


# ------------------------------
# choose a region of interest"

au = 1  # cfa index
aulab ="cfanorth"

au = 2  # cfa index
aulab ="cfasouth"

au = 3  # cfa index
aulab ="cfa4x"

yrs = 1999:2021
nP = 5  # no years to project into future 

dt = 1/12 # time resolution of solutions
no_digits = 3  # time floating point rounding 

eps = 1.0e-9 # floating point value sufficient to assume 0 valued


# ------------------------------
# dynamical model
include( "size_structured_dde!.jl" )
  if false
    include( joinpath( project_directory, "size_structured_dde_test.jl" )) # test dynamical model with generic parameters, can skip
  end


# ------------------------------
# turing model dde model
include( "size_structured_dde_turing.jl" )


# -------------------------
#  prepare data for diffeq/turing model and set default parameters
include( "size_structured_dde_turing_data.jl" )
  if false
    ## test, can ignore
    prob = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=lags )
    msol2 =  solve( prob,  solver, callback=cb, saveat=dt, isoutofdomain=(y,p,t)->any(x->x<0,y) )# to force positive
    plot( msol2, ; legend=false, xlim=(1999,2021), label="dde, with hsa, no fishing" )
  end
 

 
# ---------------
# run model estimations / overrides
 
Turing.setprogress!(false);
# Turing.setrdcache(true)

Turing.setadbackend(:forwarddiff)  # only AD that works right now
# rev diff having issues 
# Turing.setadbackend(:zygote) # 6.2 hrs, 
# Turing.setadbackend(:forwarddiff)  # CFA 4X: 6.2 hrs 
# Turing.setadbackend(:tracker)  # 5.7 hrs, but some stability issues?
# Turing.setadbackend(:reversediff)  # 5.8 hrs 

solver = MethodOfSteps(Tsit5())  
# solver = MethodOfSteps(Rodas5())  

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

 
prob = DDEProblem( size_structured_dde!, u0, h, tspan, p, constant_lags=lags )
fmod = size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nP, solver  )
 
  if false
    # for testing and timings
    # include( "size_structured_dde_turing.jl" )
    n_samples = 3
    n_adapts = 3
    n_chains = 1
    # sampler = Turing.MH()
    # sampler = Turing.HMC(0.05,10)
    sampler = Turing.NUTS(n_adapts, 0.65)
    # sampler = DynamicNUTS()
   
 
    res  =  sample( fmod, sampler, n_samples  )
 
  end

# production  
n_samples = 1000  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs
n_adapts = 1000
n_chains = Threads.nthreads()
# sampler = Turing.HMC(0.05,10)
sampler = Turing.NUTS(n_adapts, 0.65 ) # ; max_depth=10, init_ϵ=0.00625)  ;# stepsize based upon previous experience

res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )
# if on windows and threads are not working, use single processor mode:
# res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)


# ------------------------------
# save results as a hdf5

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


density(res[:"b[1]"])
density(res[:"b[2]"])
density(res[:"K[1]"])
density(res[:"v[1]"])



# -------------------------
# plot timeseries of mean fields
include( "size_structured_dde_turing_plot.jl" )

plot(0)
plots_sim = size_structured_dde_turing_plot( selection="S K predictions predictionmeans", si=[1], mw=mw ) 
# gui(plots_sim)
savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".pdf") )
savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".svg") )
savefig(plots_sim, string("size_structured_dde_turing_plots_sim_", aulab, ".png") )


plot(0)
plots_fishing = size_structured_dde_turing_plot( selection="K withfishing withoutfishing predictionmeans", si=[1], mw=mw)
# gui(plots_fishing)
savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".pdf") ) 
savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".svg") ) 
savefig(plots_fishing, string("size_structured_dde_turing_plots_fishing_", aulab, ".png") ) 



# ------------------------------
# misc computed quantities
  
pm = ( b=b, K=K, d=d, v=v, tau=tau, hsa=hsa ) 

@model function compute_derived( S1, S2, S3, S4, S5, S6, kmu, tspan, prob; nT=length(S1), M=5, er=0.2, ::Type{T}=Float64 ) where {T}  
  
  # params need to be named  .. return only that which is specified by "return()", below
  # deterministic computations: do from similations:
 
  FM = zeros(nT+M)
  BM = zeros(nT+M)
  # CA = zeros(nT+M)

  # CA[1:nT] = removed ./ K
  # CA[(nT+1):(M+nT)] = er .* bm[(nT):(M+nT-1)]
  # CA = 1.0 .- CA / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r * exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
  
  return( test=r+1, )
end

fmod2 = compute_derived(S1, S2, S3, S4, S5, S6, kmu, tspan, prob )

gq = generated_quantities(fmod2, pm)

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
 
     # Run sampler, collect results. @doc(Variational.ADVI) : 
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

     # Run sampler, collect results. @doc(Variational.ADVI) : 
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
    # Gibbs: component sampler 1, component sampler 2, ...
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


## --- rest are short snippets and notes



function prediction(x::Matrix, res, threshold)

  ## INCOMPLETE

    # Pull the means from each parameter's sampled values in the res.
    intercept = mean(res[:intercept])
    student = mean(res[:student])
    balance = mean(res[:balance])
    income = mean(res[:income])

    # Retrieve the number of rows.
    n, _ = size(x)

    # Generate a vector to store our predictions.
    v = Vector{Float64}(undef, n)

    # Calculate the logistic function for each element in the test set.
    for i in 1:n 
        num = logistic(
            intercept .+ student * x[i, 1] + balance * x[i, 2] + income * x[i, 3]
        )
        if num >= threshold
            v[i] = 1
        else
            v[i] = 0
        end
    end
    return v
end;

# Set the prediction threshold.
threshold = 0.07

# Make the predictions.
predictions = prediction(test, res, threshold)

# Calculate MSE for our test set.
loss = sum((predictions - test_label) .^ 2) / length(test_label)


    

 

#R code that needs to be ported
# frequency density of key parameters
fishery_model( DS="plot", vname="K", res=res )
fishery_model( DS="plot", vname="r", res=res )
fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 2.5))
fishery_model( DS="plot", vname="qc", res=res, xrange=c(-1, 1))
fishery_model( DS="plot", vname="FMSY", res=res  )

# timeseries
fishery_model( DS="plot", type="timeseries", vname="biomass", res=res  )
fishery_model( DS="plot", type="timeseries", vname="fishingmortality", res=res)

# Harvest control rules
fishery_model( DS="plot", type="hcr", vname="default", res=res  )

# Summary table of mean values for inclusion in document

( qs = apply(  res$mcmc$K[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )  # carrying capactiy

( qs = apply(  res$mcmc$FMSY[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) ) # FMSY


biomass = as.data.table( fit$summary("B") )
np = year.assessment+c(1:p$fishery_model$fmdata$M)
biomass$yr = rep( c(p$yrs, np ), 3)
nt = p$fishery_model$fmdata$N +p$fishery_model$fmdata$M
biomass$region = c( rep("cfanorth", nt), rep("cfasouth", nt), rep("cfa4x", nt) )
(biomass)

NN = res$p$fishery_model$fmdata$N

# densities of biomass estimates for the year.assessment
( qs = apply(  res$mcmc$B[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of biomass estimates for the previous year
( qs = apply(  res$mcmc$B[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of F in assessment year
( qs = apply(  res$mcmc$F[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,NN,], 2, mean ) )

# densities of F in previous year
( qs = apply(  res$mcmc$F[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,NN-1,], 2, mean ) )
 
 
 
 # computed quants
 
using DynamicPPL, Distributions
 @model function demo(xs)
             s ~ InverseGamma(2, 3)
             m_shifted ~ Normal(10, √s)
             m = m_shifted - 10
             for i in eachindex(xs)
                 xs[i] ~ Normal(m, √s)
             end
             return (m, )
         end
  demo (generic function with 2 methods)
  
 model = demo(randn(10));
  
 
parameters = (; s = 1.0, m_shifted=10);
  
gq = generated_quantities(model, parameters)
gq = generated_quantities(model, values(parameters), keys(parameters))


m = vec(getindex.( gq, 1))
density!(m, lab="generated quantity (VI)")
vline!([0], lab="true value")

