
# ------------------------------
# DDE size structured

# Delay Differentail EQs size structured

# https://diffeq.sciml.ai/stable/tutorials/dde_example/

# NOTE::: require 03.snowcrab_carstm.r to be completed 


pkgs = [ 
  "Revise", "MKL", "StatsBase", "Distributions", "LinearAlgebra",  "Interpolations", 
  "Plots", "StatsPlots", "MultivariateStats", "RData",
  "Turing",  "ModelingToolkit", "DifferentialEquations",  
  "StaticArrays", "LazyArrays", 
  "ForwardDiff"
  # "DiffResults", "Memoization", "DynamicPPL", "AbstractPPL", "AdvancedHMC", "MCMCChains", "SciMLSensitivity",
   #"Tracker" #, "ReverseDiff", "Zygote", "ForwardDiff", "Diffractor", "Memoization",
]
  
for pk in pkgs; @eval using $(Symbol(pk)); end

# Pkg.add( pkgs ) # add required packages

# add Turing@v0.21.9   # 21.10 error?




# Turing.setprogress!(false);
# Turing.setrdcache(true)

Turing.setadbackend(:forwarddiff)  # only AD that works right now
# rev diff having issues 
# Turing.setadbackend(:zygote) # 6.2 hrs, 
# Turing.setadbackend(:forwarddiff)  # CFA 4X: 6.2 hrs 
# Turing.setadbackend(:tracker)  # 5.7 hrs, but some stability issues?
# Turing.setadbackend(:reversediff)  # 5.8 hrs 
 
 
# include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_ode.jl")



# ------------------------------
# Part 1 -- construct basic data and parameter list defining the main characteristics of the study
 
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
        fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics" )
   
    }

    # now back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y  
    @rget Kmu 
    @rget removals 
    @rget ty
    
    # mechanism to run the rest if self contained
    include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_ode.jl")

else
    # using  RData
    
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
    o = load( fndat, convert=true)
    Y = o["Y"]
    Kmu = o["Kmu"]
    removals = o["L"]

end

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



# ------------------------------

function size_structured!( du, u, h, p, t)
  uu = max.( u, 0.0 )
  (b, K, d, v, tau, hsa)  = p
  tr21 = v[1] * max(0.0, h(p, t-1)[2])   # transition 2 -> 1   
  tr32 = v[2] * max(0.0, h(p, t-1)[3])   # transitiom 3 -> 2
  tr43 = v[3] * max(0.0, h(p, t-1)[4])   # transitiom 4 -> 3
  tr54 = v[4] * max(0.0, h(p, t-1)[5])   # transitiom 5 -> 4
  FP  = max(0.0, h(p, t-8)[6]   )      # no mature fem 8  yrs ago
  du[1] = tr21             - d[1] * uu[1] * (uu[1] / (K[1]*hsa(t,1)) )  # second order mortality       
  du[2] = tr32      - tr21 - d[2] * uu[2] * (uu[2] / (K[2]*hsa(t,2)) )  
  du[3] = tr43      - tr32 - d[3] * uu[3] * (uu[3] / (K[3]*hsa(t,3)) ) 
  du[4] = tr54      - tr43 - d[4] * uu[4] * (uu[4] / (K[4]*hsa(t,4)) ) 
  du[5] = b[1] * FP - tr54 - d[5] * uu[5] * (uu[5] / (K[5]*hsa(t,5)) )  
  du[6] = b[2] * FP        - d[6] * uu[6] * (uu[6] / (K[6]*hsa(t,6)) )   # fem mat simple logistic with lag tau and density dep on present numbers
end




if false
  # testing DDE version -- initial values for size_structured! 
  ks = 100
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* ks
  b=[ 1.0, 0.8 ]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
  d=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
  v=[0.9, 0.9, 0.9, 0.9];  
  tau=1.0 
  survey_time = 1999:2021
  forcing_time = survey_time
  external_forcing = ones(length(forcing_time),6)  # turns it off
  # external_forcing = rand(length(forcing_time),6) # random
  
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  p = ( b, K, d, v, tau, hsa )   
  tspan = (0.0, 100.0)
  nS = 6 # n components
  # history function 0.5 default
  # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5 * ks
  tau = 1  # delay
  lags = [tau]
  solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
  prob = DDEProblem( size_structured! , u0, h, tspan, p; constant_lags=lags )
  res =  solve( prob,  solver, saveat=0.1 )
  plot!( res ; legend=true)
  
  # ---------------
  #test runs
  ks = 1.0e10
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* ks
  b=[1.0, 0.8]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks; 
  d=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
  v=[0.9, 0.9, 0.9, 0.9];  
  tau=1.0 

  tspan = (1998, 2022)
  dt = 0.1

  survey_time =1999:2021

  forcing_time = survey_time
  
  external_forcing = ones(length(forcing_time),6)  # turns it off
  efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
  p = ( b, K, d, v, tau, hsa )   

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  v = 1; plot( msol.t, reduce(hcat, msol.u)'[:,v], color=[1 1] , alpha=0.75, lw=5 ) 

  plot!( msol, label="dde, no hsa, no fishing", legend=:left )
   

  external_forcing = rand(length(forcing_time),6) # random
  efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
  p = ( b, K, d, v, tau, hsa )   

  prob2 = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  # prob2 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )

  msol2 =  solve( prob2,  solver, saveat=dt )
  v = 1; plot!( msol2.t, reduce(hcat, msol2.u)'[:,v], color=[3 3] , alpha=0.75, lw=5 ) 

  plot!( msol2, label="dde, with hsa, no fishing" )

  plot!(; legend=true, xlim=(1999,2021) )
 
end




# -------------------------
# other parameters

au = 1  # cfa index
aulab ="cfanorth"

au = 2  # cfa index
aulab ="cfasouth"

au = 3  # cfa index
aulab ="cfa4x"

eps = 1.0e-9

# convert to number .. 0.56 is ave mean weight
kmu = Kmu[au] * 1000 *1000 / 0.56 * 0.9
# kmu = 1.25 * 1000 *1000 / 0.56

# "survey index"
statevars = [
  Symbol("$aulab","_M0"),
  Symbol("$aulab","_M1"),
  Symbol("$aulab","_M2"),
  Symbol("$aulab","_M3"),
  Symbol("$aulab","_M4"),
  Symbol("$aulab","_f_mat")
]
S = Matrix(Y[:, statevars ])

(nT, nS) = size(S)

dt = 0.1
yrs = 1999:2021

# spin up time of 10 years prior to start of dymamics and project 5 years into the future
tspan = (minimum(yrs)-10.0, maximum(yrs)+5.0)  

survey_time = Y[:,:yrs]   # time of observations for survey

# this only adds habitat space  ... predation is also a useful one .. 
# speed is the issue 
forcing_time = survey_time

#  sa to fraction

external_forcing =  reshape( [
   Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
   Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
   Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
   Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
   Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
   Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] ) 
  ], nT, nS )


efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
hsa = Interpolations.scale(efc, yrs, 1:nS )
 
fish_time = removals[:,:ts]
removed = removals[:,Symbol("$aulab")]

function affect_fishing!(integrator)
  i = findall(t -> t==integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
end

cb =  PresetTimeCallback( fish_time, affect_fishing! )

# history function 0.5 default
# h(p,t) = ones( nS ) .* 0.5  #values of u before t0
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* kmu
 
tau = 1  # delay
lags = [tau]

solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()


# these are dummy initial values .. just to get things started

u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* kmu
b=[1.0, 0.8]
K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu; 
d=[0.2, 0.3, 0.4, 0.5, 0.5, 0.5];
v=[0.8, 1.0, 1.0, 1.0];  
tau=1.0; 

p = ( b, K, d, v, tau, hsa )   

plot(0)

prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
msol2 =  solve( prob,  solver, saveat=dt )
plot!( msol2, label="dde, with hsa, no fishing" )

# plot!(; legend=false, xlim=(1999,2021) )


# ---------------


@model function fishery_model_turing_dde( S, kmu, tspan, prob )
    
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    nT, nS = size(S);

    K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/5.0, kmu*5.0), nS )  
    
    # consider more diffuse Cauchy prior for k .. slow mixing
    # K ~ filldist( truncated( Cauchy( kmu, kmu*0.1), kmu/10.0, kmu*10.0), nS )  
    
    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -1.0, 1.0), nS )  
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.1, 10.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.5, 0.1, 0.25, 0.9), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.5, 1.0), 4 ) 

    # initial conditions
    m = TArray{Float64}(nT, nS)
    for k in 1:nS 
      m[1,k] ~  TruncatedNormal( 0.9, 0.2, 0.5, 1.25 )  ; # starting b prior to first catch event permit higher than 1.0
    end 

    u0 = [ m[1,1], m[1,2], m[1,3], m[1,4], m[1,5], m[1,6]  ] .* K  # don't know why but takiing an array slice of m causes an error
 
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
    msol = solve( prob, solver, callback=cb, saveat=dt )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    for i in 2:nT
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        usk = max.( msol.u[j[1]], 1.0) ./ K 
        for k in 1:nS
          m[i,k] ~ TruncatedNormal( usk[k], bpsd, 1.0e-6, 1.25)  ; 
        end
      end
    end
  
    # observation model
    for k in 1:nS
      for i in 1:nT
        S[i,k] ~ TruncatedNormal( (m[i,k] + qc[k]) * q[k], bosd, 1.0e-6, 1.25 )   
      end
    end
    
end

 
# ---------------


prob = DDEProblem( size_structured!, u0, h, tspan, p, constant_lags=lags )
fmod = fishery_model_turing_dde( S, kmu, tspan, prob )
# fmod = fishery_model_turing_incremental_dde( S, kmu,  tspan, prob )


# testing
n_samples = 3
n_adapts = 3
n_chains = 1
sampler = Turing.MH()
sampler = Turing.HMC(0.05,10)
# sampler = Turing.NUTS(n_adapts, 0.65)
res  =  sample( fmod, sampler, n_samples  )


# production .. ~ 5 hrs 
n_samples = 500
n_adapts = 500
n_chains = Threads.nthreads()
sampler = Turing.NUTS(n_adapts, 0.65)


res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )
# if on windows o threads not working:
# res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)


# save as a hdf5
using JLD2  # using HDF5
fn = string("/home/jae/julia/snowcrab/data_size_struct_dde", "_", aulab, ".hdf5" )
@save fn res
@load fn res

# in R:  h5read( paste("/home/jae/julia/snowcrab/data_size_struct_dde", "_", aulab, ".hdf5"), "res")

# native julaia dat format
fn = string("/home/jae/julia/snowcrab/data_size_struct_dde", "_", aulab, ".jd" )
write(fn, res)
read(fn)
 
 
show(stdout, "text/plain", summarize(res))

density(res[:"b[1]"])
density(res[:"b[2]"])
density(res[:"K[1]"])
density(res[:"v[1]"])

# rng = MersenneTwister(26)
# resp = predict(rng, textmodel_marginal_pred(data), res)
# @df resp ecdfplot(:"b[1]"; label="birth rate 1")


# mean field dynamics:
u0 = [ 
  mean( res[:,"K[1]",:] ),
  mean( res[:,"K[2]",:] ),
  mean( res[:,"K[3]",:] ),
  mean( res[:,"K[4]",:] ),
  mean( res[:,"K[5]",:] ),
  mean( res[:,"K[6]",:] )
] .*  [ 
  mean( res[:,"m[1,1]",:] ),
  mean( res[:,"m[1,2]",:] ),
  mean( res[:,"m[1,3]",:] ),
  mean( res[:,"m[1,4]",:] ),
  mean( res[:,"m[1,5]",:] ),
  mean( res[:,"m[1,6]",:] )
]

b = [ mean( res[:,"b[1]",:] ), mean( res[:,"b[2]",:] ) ]
K = [ mean( res[:,"K[1]",:] ), mean( res[:,"K[2]",:] ), 
      mean( res[:,"K[2]",:] ), mean( res[:,"K[2]",:] ),
      mean( res[:,"K[2]",:] ), mean( res[:,"K[2]",:] ) ]  ; 
d = [ mean( res[:,"d[1]",:] ), mean( res[:,"d[2]",:] ), 
      mean( res[:,"d[3]",:] ), mean( res[:,"d[4]",:] ),
      mean( res[:,"d[5]",:] ), mean( res[:,"d[6]",:] ) ]   
v = [ mean( res[:,"v[1]",:] ), mean( res[:,"v[2]",:] ), 
      mean( res[:,"v[3]",:] ), mean( res[:,"v[4]",:] ) ]

pm = ( b, K, d, v, tau, hsa ) 

msol = solve( remake( prob, u0=u0, h=h, tspan=tspan, p=pm ), solver, callback=cb, saveat=dt )  
plot!(msol, label="dde-mean-field-fishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[1 v] , alpha=0.5, lw=5 ) 
plot!(; legend=false, xlim=(1997,2023) )

prob2 = DDEProblem( size_structured!, u0, h, tspan, pm, saveat=dt )
msol = solve( prob2,  solver, saveat=dt ) #  effective nullify callbacks
plot!(msol, label="dde-mean-field-nofishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[5 5] , alpha=0.75, lw=5 ) 


# back transform S to normal scale 
j = 1  # fishable component
yhat = ( S[:,j] .* mean(res[:,Symbol("q[$j]"),:]) .- mean(res[:,Symbol("qc[$j]"),:] ) ) .* mean(res[:,Symbol("K[$j]"),:] ) 
scatter!(survey_time, yhat   ; color=[1 2])
plot!(survey_time, yhat  ; color=[j j])
plot!(; legend=false, xlim=(1997,2023) )


# sample and plot posterior K
j = 1  # state variable index
for u in 1:length(res)  
  w = res[u,Symbol("K[$j]"),1]
  hline!([w];  alpha=0.1, color=[j j])
end

plot!(; legend=false, xlim=(1997,2023) )


# sample and plot posterior means from model (posterior post-fishery abundance)
j = 1  # state variable index
#j = 6
w = zeros(nT)
for u in 1:length(res)  
  for i in 1:nT
    w[i] = res[u,Symbol("K[$j]"),1] * res[u,Symbol("m[$i,$j]"),1]
  end
  plot!(survey_time, w  ;  alpha=0.1, color=[j j])
end

plot!(; legend=false, xlim=(1997,2023) )

# mean post-fishery abundance
u = zeros(nT)
v = zeros(nT)
for  i in 1:nT
  u[i] = mean( res[:,Symbol("m[$i,$j]"),:] .* res[:,Symbol("K[$j]"),:] ) 
  v[i] = std(  res[:,Symbol("m[$i,$j]"),:] .* res[:,Symbol("K[$j]"),:] ) 
end
scatter!(survey_time, u  ; color=[3 2])

  
# misc computed quantities

# params need to be named  .. return only that which is specified by "return()", below
pm = ( b=b, K=K, d=d, v=v, tau=tau, hsa=hsa ) 

@model function fm_test( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, nT=length(S1), ::Type{T}=Float64 ) where {T}  
  
  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(nT+M)
  B = zeros(nT+M)
  C = zeros(nT+M)

  C[1:nT] = removed ./ K
  C[(nT+1):(M+nT)] = er .* bm[(nT):(M+nT-1)]
  C = 1.0 .- C / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = bm .* K  
 
  # recaled estimates
  B[1:nT] = bm[1:nT] *. K - L[1:nT] ;
  B[(nT+1):(M+nT)] = (bm[(nT+1):(M+nT)] - C[(nT):(M+nT-1)]) *. K ;
 
  return( test=r+1, )
end

fmod2 = fm_test(S1, S2, S3, S4, S5, S6, kmu, tspan, prob )

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

msol =  solve( prob,  solver, callback=cb, saveat=dt)  

function predict_rd() # Our 1-layer "neural network"
  solve(prob, solver,p=p,saveat=dt)[1] # override with new parameters
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


