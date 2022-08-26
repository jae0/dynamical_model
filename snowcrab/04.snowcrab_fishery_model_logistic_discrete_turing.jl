
# ------------------------------
# Discrete

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


for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages

 
Turing.setprogress!(false);
# Turing.setrdcache(true)

Turing.setadbackend(:forwarddiff)   
  # Turing.setadbackend(:forwarddiff)  #  
  # Turing.setadbackend(:zygote) #  
  # Turing.setadbackend(:tracker)  #  
  # Turing.setadbackend(:reversediff)  #  
  


# ------------------------------
# choose a region of interest"

au = 1  # cfa index
aulab ="cfanorth"

au = 2  # cfa index
aulab ="cfasouth"

au = 3  # cfa index
aulab ="cfa4x"

yrs = 1999:2021  # <<<<<<<<-- change
nP = 5  # number of predictions into future (with no fishing)

dt = 1  # time resolution of solutions
no_digits = 3  # time floating point rounding 

eps = 1.0e-9 # floating point value sufficient to assume 0 valued

 
# ------------------------------
# turing model -- dynamics embedded
include( "logistic_discrete_turing.jl" )


# -------------------------
#  prepare data for turing model and set default parameters
include( "logistic_discrete_turing_data.jl" )
 

# ---------------
# run model estimations / overrides
Turing.setprogress!(false);
 


# ------------------------------
# initial values for Logistic!
# almost working but solution decay to negative numbers though it is supposed to be bounded ..  

# choose one:    
model_variation = "Model_1"
fmod = logistic_discrete_turing_basic( S, kmu, nT, nP, removed  )  # q only

model_variation = "Model_2"
fmod = logistic_discrete_turing( S, kmu, nT, nP, removed  ) # q, qc

  if false 
    # Prior predictive check: 
    prior_res = sample(   fmod, Prior(), 100, nwarmup = 100, nchains = 3 );

    missing_data = Vector{Missing}(missing, nT)
    mod_preds = fmod( missing_data, kmu, nT, removed )
    prior_check = predict( mod_preds, prior_res )
    summarystats( prior_check )
    plot(prior_check)
  end

# Posterior sampling
  if false
    # for testing and timings
    # include( "logistic_discrete_turing.jl" )
    n_samples = 3
    n_adapts = 3
    n_chains = 1
    # sampler = Turing.MH()
    # sampler = Turing.HMC(0.05,10)
    sampler = Turing.NUTS(n_adapts, 0.65; max_depth=10, init_ϵ=0.05)
    # sampler = DynamicNUTS()

    res  =  sample( fmod, sampler, n_samples  )
  end

# production  
n_samples = 2000
n_adapts = 1000
n_chains = Threads.nthreads()
# sampler = Turing.HMC(0.05,10)
sampler = Turing.NUTS(n_adapts, 0.95; max_depth=12, init_ϵ=0.05)  ;# stepsize based upon previous experience

res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )  # < 1 min
# if on windows and threads are not working, use single processor mode:
# res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)



# ------------------------------
# save results as a hdf5

fn = joinpath( project_directory, string("logistic_discrete_turing_data", "_", model_variation, "_", aulab, ".hdf5" ) )
@save fn res
# @load fn res
# can read back in R as:  
# h5read( paste("/home/jae/julia/snowcrab/logistic_discrete_turing_data","_", model_variation, "_", aulab, ".hdf5"), "res")

# fnk = joinpath( project_directory, string("logistic_discrete_turing_data", "_", model_variation,"_", aulab, "_K_", ".hdf5" ) )
# tt =res[:,Symbol("K[1]"),:]
# @save fnk tt

# ------------------------------
# process outputs
# display all estimates
show(stdout, "text/plain", summarize(res))


density(res[:"q"])
density(res[:"qc"])
density(res[:"K"])
density(res[:"r"])


  show(stdout, "text/plain", summarize(res))

  histogram(res[:"r"])
  histogram(res[:"K"])
  histogram(res[:"q"])
  histogram(res[:"qc"])
  histogram(res[:"bosd"])
  histogram(res[:"bpsd"])

 
# -------------------------
# plot timeseries of mean fields
include( "logistic_discrete_turing_plot.jl" )

plot(0)
plots_sim = logistic_discrete_turing_plot( selection="S K predictions predictionmeans"  ) 
# gui(plots_sim)
savefig(plots_sim, string("logistic_discrete_turing_plots_sim", "_", model_variation, "_", aulab, ".pdf") )
savefig(plots_sim, string("logistic_discrete_turing_plots_sim", "_", model_variation, "_", aulab, ".svg") )
savefig(plots_sim, string("logistic_discrete_turing_plots_sim", "_", model_variation, "_", aulab, ".png") )


plot(0)
plots_fishing = logistic_discrete_turing_plot( selection="K predictionmeans" )
# gui(plots_fishing)
savefig(plots_fishing, string("logistic_discrete_turing_plots_fishing", "_", model_variation, "_", aulab, ".pdf") ) 
savefig(plots_fishing, string("logistic_discrete_turing_plots_fishing", "_", model_variation, "_", aulab, ".svg") ) 
savefig(plots_fishing, string("logistic_discrete_turing_plots_fishing", "_", model_variation, "_", aulab, ".png") ) 


    
  using MCMCChains, AdvancedHMC

  summarystats(res)

  density( res[:,:r,:])
  density( res[:,:K,:])
  density( res[:,:q,:])
  density( res[:,:qc,:])
  

  group( res, :m0)  ##== res(:, 5:20, :) == res[[:m0]]
  
  corner(res)

  plot( res[[:m0]] )

  plot(
    traceplot(res),
    meanplot(res),
    density(res),
    histogram(res),
    mixeddensity(res),
    autocorplot(res),
    dpi=300, size=(840,600)
  )

  plot(res, seriestype=(:meanplot, :autocorplot), dpi=300)


  plot(res[:,[Symbol("m0[$i]") for i in 1:10],:])

  using ArviZ
  using PyPlot



# ArviZ ships with style sheets!
ArviZ.use_style("arviz-darkgrid")
  

plot_autocorr(res; var_names=["r", "K"]);
gcf()

idata = from_mcmcchains( res; library="Turing" )

Plots.plot( survey_time , summarystats(idata.posterior; var_names=["m0"]).mean )
Plots.plot!( survey_time , S )


Plots.plot( survey_time , summarystats(idata.posterior; var_names=["ys"]).mean )



oo = Array(res, (size(res)[1], size(res)[2]) )

  kmean = mean( res[[3]].value .*  res[[2]].value)

  pm = [mean( res[[1]].value ), mean( res[[2]].value ) ]

  odeprob = ODEProblem( Logistic!, [kmean], tspan, pm )
  bm = solve( odeprob, Rosenbrock23(); saveat=0.1, callback=cb, reltol=1e-15,abstol=1e-15 )
  plot!(bm, label="ode-fishing")

  prob = remake( prob, u0=[kmean], p=pm)
  bm = solve(prob, SSAStepper(),  saveat=0.1, callback=cb  ) #
  plot!(bm, lw=2, label="jump-fishing" )
 
  plot(; legend=false)
  posterior_samples = sample(res[[1, 2, 3]], 50; replace=false)
  for p in eachrow(Array(posterior_samples))
      bmpost = solve(prob, Tsit5(); p=pm, saveat=0.1)
      plot!(bmpost; alpha=0.1, color="#BBBBBB")
  end

  # Plot simulation and noisy observations.
  bm =  solve( prob, AutoTsit5(Rosenbrock23()); saveat=0.1, callback=cb, reltol=1e-16, abstol=1e-16 ) 
  plot!(bm; color=[1 2], linewidth=1)
  scatter!(bm.t, S'; color=[1 2])


end



  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(sN+M)
  B = zeros(nT+M)
  C = zeros(nT+M)

  C[1:nT] = removals ./ K
  C[(nT+1):(M+nT)] = er .* bm[(nT):(M+nT-1)]
  C = 1.0 .- C / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = bm .* K  
 



bm = tzeros(Real, nT+M, U)
  
bm[1] ~ TruncatedNormal( b0, bpsd, eps, 1.0)  ;
for i in 2:nT 
  o = r * ( 1.0 - bm[i-1] ) ; 
  bm[i] ~ TruncatedNormal(   bm[i-1] * ( 1.0 + o ) - L[i-1]/K, bpsd, eps, 1.0)  ;
end
for i in (nT+1):(M+nT) 
  o = r * ( 1.0 - bm[i-1] ) ; 
  bm[i] ~ TruncatedNormal(   bm[i-1] * ( 1.0 + o ) - er*bm[(i-1)], bpsd, eps, 1.0)  ;
end

# -------------------
# biomass observation model
# cfanorth(1) and cfasouth(2)
#   This is slightly complicated because a fall / spring survey correction is required:
#   B represents the total fishable biomass available in fishing year y
#     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
#     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
# spring surveys from 1998 to 2003
#   this is conceptualized in the following time line:
#     '|' == start/end of each new fishing year
#     Sf = Survey in fall
#     Ss = Survey in spring
#     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
# Cfa 4X -- fall/winter fishery
#    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


# north and south
if j in 1:2 {
  # spring surveys
  ys = ( Y[1, j] / q ) +  qc
  ys ~ TruncatedNormal( bm[1,j] - L[1,j]/K , bosd, eps, 1.0) ;
  for i in 2:(ty-1) 
    ys = ( Y[i, j] / q ) +  qc
    ys  ~ TruncatedNormal( bm[i,j] - L[i-1,j]/K , bosd, eps, 1.0)  ;
  end
  #  transition year (ty)
  ys = ( Y[ty,j] / q ) +  qc
  ys  ~ TruncatedNormal(  bm[ty,j]  - (L[ty-1,j]/K  + L[ty,j]/K ) / 2.0  , bosd, eps, 1.0)  ; #NENS and SENS
  # fall surveys
  for j in 1:U 
    for i in (ty+1):nT 
      ys = ( Y[i,j] / q ) +  qc
      ys ~ TruncatedNormal(  bm[i,j] - L[i,j]/K, bosd, eps, 1.0)  ; #   fall surveys
    end
  end

end

# cfa4X
if j ==3
  # spring surveys
  for i in 1:(ty-1)  
    ys = ( Y[i, 3] / q[3] ) +  qc[3]
    ys  ~ TruncatedNormal( bm[i,3] - L[i,3]/K[3], bosd[3], eps, 1.0)  ;
  end
  #  transition year (ty)
  ys = ( Y[ty,3] / q[3] ) +  qc[3]
  ys  ~ TruncatedNormal(  bm[ty,3]  - L[ty,3]/K[3] , bosd[3], eps, 1.0)  ; #SENS
  # fall surveys
  for j in 1:U 
    for i in (ty+1):nT 
      ys = ( Y[i,j] / q ) +  qc
      ys ~ TruncatedNormal(  bm[i,j] - L[i,j]/K, bosd, eps, 1.0)  ; #   fall surveys
    end
  end

end  

# deterministic computations: 
F = zeros(sN+M)
B = zeros(nT+M)
C = zeros(nT+M)

C[1:nT] = L[1:nT] ./ K
C[(nT+1):(M+nT)] = er .* bm[(nT):(M+nT-1)]
C = 1.0 -. C / bm

F =  -log( max.(C, eps) )  ;

   
# -------------------
# parameter estimates for output
MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
BMSY   = exp(K)/2 ; # biomass at MSY
FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY


# recaled estimates

B[1:nT] = bm[1:nT] *. K - L[1:nT] ;
B[(nT+1):(M+nT)] = (bm[(nT+1):(M+nT)] - C[(nT):(M+nT-1)]) *. K ;


 
