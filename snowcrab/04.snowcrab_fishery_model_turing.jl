# ----------------
# example template calling julia 
project_directory = string(expanduser("~/projects/dynamical_model/"), "snowcrab")
# project_directory = @__DIR__() #  same folder as the file

push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(project_directory)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  
  "Plots", "StatsPlots", "MultivariateStats"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages


# Part 1 -- construct basic parameter list defining the main characteristics of the study

# NOTE::: require 03.snowcrab_carstm.r to be completed 
# (i.e.,spatiotemporal model and assimilate_numbers_and_size to have been completed 

get_data_with_RCall = false

if get_data_with_RCall

    using RCall

    # typing <$> in Julia's  command prompt starts an R session.  
    
    $
      
    {
        # this is R-code that creates local RData file with required data
        source( file.path( code_root, "bio_startup.R" )  )
        require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics" ) 
        fishery_model_data_inputs( year.assessment=2021, type="numerical_dynamics" )
        fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics" )
        # type <backspace> to escape back to julia
    }

    # now back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y  
    @rget Kmu 
    @rget Ksd 
    @rget removals 
    @rget ty
    
    # mechanism to run the rest if self contained
    include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_ode.jl")

else

    using CodecBzip2, CodecXz, RData  
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_biomass.RData"
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number.RData"
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
    o = load( fndat, convert=true)
    Y = o["Y"]
    Ksd = o["Ksd"]
    Kmu = o["Kmu"]
    removals = o["L"]

end



Turing.setprogress!(false);
Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
# Turing.setadbackend(:tracker)
 
  
 
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

