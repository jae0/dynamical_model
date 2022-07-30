
# ------------------------------
# Discrete

# NOTE::: require 03.snowcrab_carstm.r to be completed 



project_directory = string(expanduser("~/projects/dynamical_model/"), "snowcrab")
# project_directory = @__DIR__() #  same folder as the file

push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(project_directory)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions", "DynamicPPL",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  "Interpolations", 
  "Plots", "StatsPlots", "MultivariateStats", "RData"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages





# ------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
 
 
get_data_with_RCall = false

if get_data_with_RCall

      
    {
        # this is R-code that creates local RData file with required data
        # run in R externally or from within julia or .. 
        { # from within julia
          using RCall
          # type $ in Julia's command prompt starts an R session.  
          # .. run below
          # type <backspace> to escape back to julia
        }
    
        source( file.path( code_root, "bio_startup.R" )  )
        require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics", for_julia=TRUE )
   
    }

    # now back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y  
    @rget Kmu 
    @rget removals 
    @rget ty
    
    # mechanism to run the rest if self contained
        fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics" )
    }

    # now back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y  
    @rget Kmu 
    @rget Ksd 
    @rget removals 
    @rget ty
    
    # mechanism  to run the rest if self contained
    include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_discrete.jl")

else

    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_biomass.RData"
    o = load( fndat, convert=true)
    Y = o["Y"]
    Ksd = o["Ksd"]
    Kmu = o["Kmu"]
    removals = o["L"]

end


# ------------------------------

  # initial values for Logistic!
  # almost working but solution decay to negative numbers though it is supposed to be bounded ..  
    
  cfanorth =  1 # column index
  cfasouth =  2 # column index
  cfa4x =  3 # column index

  au = cfanorth  # cfa index
  aulab ="cfanorth"
  eps = 1.0e-9
  
  # convert to number .. 0.56 is ave mean weight of fb
  kmu = Kmu[au]  
  ksd = Ksd[au]

  
  # "survey index"
  S1 = Y[:,Symbol("$aulab"  )]

  nS = 1 # n components
  nT = length(S1)
  dt = 1
  yrs = 1999:2021
  
  survey_time = Y[:,:yrs]   # time of observations for survey
 
  fish_time = removals[:,:ts]
  removed = removals[:,Symbol("$aulab")]


  

  @model function fishery_model_turing_discrete( S1, kmu, removed, ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    eps = 1.0e-9
    nT=length(S1)

    K ~ TruncatedNormal( kmu, kmu*0.1, kmu/2.5, kmu*2.5)  
    r ~  TruncatedNormal( 1.0, 0.1, 0.25, 1.75)   # (mu, sd)

    bosd ~ TruncatedNormal( 0.1, 0.05, eps, 0.25)  ;  # slightly informative .. center of mass between (0,1)
    bpsd ~ TruncatedNormal( 0.1, 0.05, eps, 0.25)  ;  # slightly informative .. center of mass between (0,1)

    q ~ TruncatedNormal( 1.0, 0.05, 0.5, 2.0)    
    qc ~ TruncatedNormal(0.0, 0.05, -1.0, 1.0) 

    # m1 = fished (postfishery) abundance
    m1 =  Vector{T}(undef, nT)
    m1[1] ~  TruncatedNormal( 0.9, 0.05, eps, 1.25)  ; # starting b prior to first catch event

    for i in 2:nT
      m1[i] ~ TruncatedNormal( r * m1[i-1] * ( 1.0 - m1[i-1] ) - removed[i-1]/K, bpsd, eps, 1.25)  ;
    end

    ys = Vector{T}(undef, nT)
    for i in 1:nT
      S1[i] ~ TruncatedNormal( (m1[i] - qc)*q, bosd, eps, 1.25 )  ;
      ys[i] = m1[i] * K
    end
  end


  fmod = fishery_model_turing_discrete( S1, kmu, removed  )
   
  
  # Prior predictive check: 
  prior_res = sample(   fmod, Prior(), 100, nwarmup = 100, nchains = 3 );
 
  missing_data = Vector{Missing}(missing, nT)
  mod_preds = fmod( missing_data, kmu, removed )
  prior_check = predict( mod_preds, prior_res )
  summarystats( prior_check )
  plot(prior_check)

  # Posterior sampling

  # testing
  n_samples = 3
  n_adapts = 3
  n_chains = 1
  sampler = Turing.MH()
  sampler = Turing.NUTS(n_adapts, 0.65)
  res  =  sample( fmod, sampler, n_samples )


  # production
  n_samples = 2000
  n_adapts = 1000
  n_chains = 3
  sampler = Turing.NUTS(n_adapts, 0.95)


  res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )
  # if on windows o threads not working:
  # res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)

  show(stdout, "text/plain", summarize(res))

  histogram(res[:"r"])
  histogram(res[:"K"])
  histogram(res[:"q"])
  histogram(res[:"qc"])

  # sample and plot means from model

  w = zeros(nT)
  for u in 1:length(res)  
    for i in 1:nT
      w[i] = res[u,Symbol("K"),1] * res[u,Symbol("m1[$i]"),1]
    end
    plot!(survey_time, w  ;  alpha=0.05, color=[2 2])
  end
  plot!(; legend=false, xlim=(1997,2023), ylim=(0,10) )

  u = zeros(nT)
  v = zeros(nT)
  for  i in 1:nT
    u[i] = mean( res[:,Symbol("m1[$i]"),:] .* res[:,Symbol("K"),:] ) 
    v[i] = std(  res[:,Symbol("m1[$i]"),:] .* res[:,Symbol("K"),:] ) 
  end
  scatter!(survey_time, u;   alpha=0.75, color=[3 2])

    
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
Plots.plot!( survey_time , S1 )


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
  scatter!(bm.t, S1'; color=[1 2])


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


 
