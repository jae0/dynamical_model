


# ------------------------------
# DDE size structured + SSA


# Delay Differentail EQs size structured

# https://diffeq.sciml.ai/stable/tutorials/dde_example/

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
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  "Interpolations", "JumpProcesses", "DelaySSAToolkit",
  "Plots", "StatsPlots", "MultivariateStats"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages

Threads.nthreads()



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

Turing.setprogress!(false);
# Turing.setadbackend(:zygote)
Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
# Turing.setadbackend(:tracker)
 

# markovian part ... use: ode or ssa


 
# nonmarkovian .. delayssatoolkit or delaydifferential equation

function size_structured!( du, u, h, p, t)
  # S1, S2, S3, S4, S5, S6 = u
  b, K, d, v, tau, hsa  = p
  tr21 = v[1] * h(p, t-1)[2]   # transition 2 -> 1   
  tr32 = v[2] * h(p, t-1)[3]   # transitiom 3 -> 2
  tr43 = v[3] * h(p, t-1)[4]   # transitiom 4 -> 3
  tr54 = v[4] * h(p, t-1)[5]   # transitiom 5 -> 4
  f8  = h(p, t-8)[6]           # no fem 8 yrs ago
  du[1] = tr21             - (d[1] * u[1]) * (u[1]/ K[1]) * hsa(t,1)       
  du[2] = tr32      - tr21 - (d[2] * u[2]) * (u[2]/ K[2]) * hsa(t,2) 
  du[3] = tr43      - tr32 - (d[3] * u[3]) * (u[3]/ K[3]) * hsa(t,3)
  du[4] = tr54      - tr43 - (d[4] * u[4]) * (u[4]/ K[4]) * hsa(t,4)
  du[5] = b[1] * f8 - tr54 - (d[5] * u[5]) * (u[5]/ K[5]) * hsa(t,5) 
  du[6] = b[2] * f8        - (d[6] * u[6]) * (u[6]/ K[6]) * hsa(t,6)  # fem mat simple logistic with lag tau and density dep on present numbers
end




if false
  # testing DDE version -- initial values for size_structured! 
  
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 100
  b=[  3.0, 1.24]
  K=[100.0, 100.0, 100.0, 100.0, 100.0, 100.0];
  d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
  v=[0.8, 0.8, 0.8, 0.8];  
  tau=1.0 

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
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
  tau = 1  # delay
  lags = [tau]
  solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
  prob = DDEProblem( size_structured , u0, h, tspan, p; constant_lags=lags )
  res =  solve( prob,  solver, saveat=0.1 )
  plot( res )
 

  # ---------------
  #test runs

  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 1.e9
  b=[3.0, 1.24]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*1.e-9; 
  d=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
  v=[0.9, 0.9, 0.9, 0.8];  
  tau=1.0; 

  p = ( b, K, d, v, tau, hsa )   

  forcing_time = survey_time

  external_forcing = ones(length(forcing_time),6)  # turns it off
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  plot!( msol, label="dde, no hsa, no fishing", legend=:left )
   

  external_forcing = rand(length(forcing_time),6) # random
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  plot!( msol, label="dde, with hsa, no fishing", legend=:left )


end




# -------------------------
# other parameters

au = 2  # cfa index
aulab ="cfasouth"
eps = 1.0e-9

# convert to number .. 0.56 is ave mean weight of fb
kmu = Kmu[au] * 1000 *1000 / 0.56

# "survey index"
S1 = Y[:,Symbol("$aulab","_M0")]  
S2 = Y[:,Symbol("$aulab","_M1")]  
S3 = Y[:,Symbol("$aulab","_M2")]  
S4 = Y[:,Symbol("$aulab","_M3")] 
S5 = Y[:,Symbol("$aulab","_M4")]  
S6 = Y[:,Symbol("$aulab","_f_mat")]  

nS = 6 # n components
nT = length(S1)
dt = 0.1
yrs = 1999:2021
tspan = (1999.0, 2021.0)

survey_time = Y[:,:yrs]   # time of observations for survey

# this only adds habitat space  ... predation is also a useful one .. 
# speed is the issue 
forcing_time = survey_time

# invert sa to fraction

external_forcing = reshape( [
    1.0 .- Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] ) 
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
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5 *kmu
 
tau = 1  # delay
lags = [tau]

solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()


# these are dummy initial values .. just to get things started
u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* kmu
b=[1.0, 1.0]
K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*1.e-9; 
d=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
v=[0.9, 0.9, 0.9, 0.8];  
tau=1.0; 

p = ( b, K, d, v, tau, hsa )   

 
# ---------------


@model function fishery_model_turing_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, nS=6, N=length(S1), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K ~ filldist( truncated( Cauchy( kmu, kmu*0.25), kmu/10.0, kmu*10.0), nS )  

    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ filldist( truncated( Cauchy( 1.0, 1.0), 0.1, 10.0), nS )    
    qc ~ filldist( truncated( Cauchy(0.0, 0.1), -1.0, 1.0), nS )  
  
    # initial conditions
    m1 =  Vector{T}(undef, N)
    m2 =  Vector{T}(undef, N)
    m3 =  Vector{T}(undef, N)
    m4 =  Vector{T}(undef, N)
    m5 =  Vector{T}(undef, N)
    m6 =  Vector{T}(undef, N)

    m1[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m2[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m3[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m4[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m5[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m6[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event

   # birth rate from F_8
    b ~ filldist( truncated( Cauchy(1.0, 0.1), 0.5, 2.0), 2 ) 
     
    # mortality
    d ~ filldist( truncated( Cauchy(0.2, 0.1), 1.0e-9, 0.9), nS )  

    # transition rates
    v ~ filldist( truncated( Cauchy(0.9, 0.1), 1.0e-9, 1.0 ), 4 ) 
  
    # process model
    u0 = T[ m1[1], m2[1], m3[1], m4[1], m5[1], m6[1] ] .* K 
    p = ( b, K, d, v, tau, hsa )
    msol = solve( remake( prob; u0=u0, h=h, tspan=tspan, p=p ), solver, callback=cb, saveat=dt )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        usk = max.( msol.u[j[1]] ./ K, 1.0e-9 )
        m1[i] ~ TruncatedNormal( usk[1], bpsd, 1.0e-9, 1.0)  ; 
        m2[i] ~ TruncatedNormal( usk[2], bpsd, 1.0e-9, 1.0)  ; 
        m3[i] ~ TruncatedNormal( usk[3], bpsd, 1.0e-9, 1.0)  ; 
        m4[i] ~ TruncatedNormal( usk[4], bpsd, 1.0e-9, 1.0)  ; 
        m5[i] ~ TruncatedNormal( usk[5], bpsd, 1.0e-9, 1.0)  ; 
        m6[i] ~ TruncatedNormal( usk[6], bpsd, 1.0e-9, 1.0)  ; 
      end
    end
  
    # observation model
    @. S1 ~ TruncatedNormal( (m1 + qc[1]) * q[1], bosd, 0.0, 1.0 )   
    @. S2 ~ TruncatedNormal( (m2 + qc[2]) * q[2], bosd, 0.0, 1.0 )   
    @. S3 ~ TruncatedNormal( (m3 + qc[3]) * q[3], bosd, 0.0, 1.0 )   
    @. S4 ~ TruncatedNormal( (m4 + qc[4]) * q[4], bosd, 0.0, 1.0 )   
    @. S5 ~ TruncatedNormal( (m5 + qc[5]) * q[5], bosd, 0.0, 1.0 )   
    @. S6 ~ TruncatedNormal( (m6 + qc[6]) * q[6], bosd, 0.0, 1.0 )   
    
end


# ---------------

@model function fishery_model_turing_incremental_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, N=length(S1), ::Type{T}=Float64 ) where {T}  
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors
  K ~ filldist( truncated( Cauchy( kmu, kmu*0.25), kmu/10.0, kmu*10.0), nS )  

  bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ filldist( truncated( Cauchy( 1.0, 1.0), 0.1, 10.0), nS )    
  qc ~ filldist( truncated( Cauchy(0.0, 0.1), -1.0, 1.0), nS )  

  # initial conditions
  m1 =  Vector{T}(undef, N)
  m2 =  Vector{T}(undef, N)
  m3 =  Vector{T}(undef, N)
  m4 =  Vector{T}(undef, N)
  m5 =  Vector{T}(undef, N)
  m6 =  Vector{T}(undef, N)

  m1[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
  m2[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
  m3[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
  m4[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
  m5[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
  m6[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event

  # birth rate from F_8
  b ~ filldist( truncated( Cauchy(1.0, 0.1), 0.25, 3.0), 2 ) 
   
  # mortality
  d ~ filldist( truncated( Cauchy(0.2, 0.1), 1.0e-9, 0.9), nS )  

  # transition rates
  v ~ filldist( truncated( Cauchy(0.9, 0.1), 1.0e-9, 1.0 ), 4 ) 
 
  # process model
  t0 = floor(survey_time[1]) 
  
  for i in 2:N
    tsp = (t0+i-1.1, t0+i+0.1 )
    u0 = T[ m1[1], m2[1], m3[1], m4[1], m5[1], m6[1] ] .* K 
    p = ( b, K, d, v, tau, hsa )
  
    msol = solve( 
      remake( prob; u0=u0, h=h, tspan=tsp, p=p, constant_lags=lags   ), 
      solver, 
      callback=cb,  
      saveat=dt 
    )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    j = findall(t -> t==survey_time[i], msol.t)
    if length(j) > 0
      usk = max.( msol.u[j[1]] ./ K, 1.0e-9 )
      m1[i] ~ TruncatedNormal( usk[1], bpsd, 1.0e-9, 1.0)  ; 
      m2[i] ~ TruncatedNormal( usk[2], bpsd, 1.0e-9, 1.0)  ; 
      m3[i] ~ TruncatedNormal( usk[3], bpsd, 1.0e-9, 1.0)  ; 
      m4[i] ~ TruncatedNormal( usk[4], bpsd, 1.0e-9, 1.0)  ; 
      m5[i] ~ TruncatedNormal( usk[5], bpsd, 1.0e-9, 1.0)  ; 
      m6[i] ~ TruncatedNormal( usk[6], bpsd, 1.0e-9, 1.0)  ; 
    end
  end

  # observation model
  @. S1 ~ TruncatedNormal( (m1 + qc[1]) * q[1], bosd, 0.0, 1.0 )   
  @. S2 ~ TruncatedNormal( (m2 + qc[2]) * q[2], bosd, 0.0, 1.0 )   
  @. S3 ~ TruncatedNormal( (m3 + qc[3]) * q[3], bosd, 0.0, 1.0 )   
  @. S4 ~ TruncatedNormal( (m4 + qc[4]) * q[4], bosd, 0.0, 1.0 )   
  @. S5 ~ TruncatedNormal( (m5 + qc[5]) * q[5], bosd, 0.0, 1.0 )   
  @. S6 ~ TruncatedNormal( (m6 + qc[6]) * q[6], bosd, 0.0, 1.0 )   
  
end

# ---------------

prob = DDEProblem( size_structured!, u0, h, tspan, p, constant_lags=lags )
fmod = fishery_model_turing_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob )
# fmod = fishery_model_turing_incremental_dde( S1, S2, S3, S4, S5, S6, kmu,  tspan, prob )


# testing
n_samples = 3
n_adapts = 3
n_chains = 1
sampler = Turing.MH()
sampler = Turing.NUTS(n_adapts, 0.65)


# production
n_samples = 1000
n_adapts = 500
n_chains = 3
sampler = Turing.NUTS(n_adapts, 0.9)


res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )
# if on windows o threads not working:
# res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)

show(stdout, "text/plain", summarize(res))

histogram(res[:"b[1]"])
histogram(res[:"b[2]"])


t0 = floor(survey_time[1])

j = 1 # S1 

plot(; legend=false, xlim=(1997,2023) )

for u in 1:n_samples
  
  for  i in 1:nT
    u0 = [ 
      res[u,Symbol("m1[$i]"),1] * res[u,:"K[1]",1],
      res[u,Symbol("m2[$i]"),1] * res[u,:"K[2]",1],
      res[u,Symbol("m3[$i]"),1] * res[u,:"K[3]",1],
      res[u,Symbol("m4[$i]"),1] * res[u,:"K[4]",1],
      res[u,Symbol("m5[$i]"),1] * res[u,:"K[5]",1],
      res[u,Symbol("m6[$i]"),1] * res[u,:"K[6]",1]
    ]
    
    tspn = ( t0+i-1.1, t0+i+0.1 )

    b = [ res[u,:"b[1]",1], res[u,:"b[2]",1] ]
    K = [ res[u,:"K[1]",1], res[u,:"K[2]",1], res[u,:"K[3]",1], res[u,:"K[4]",1], res[u,:"K[5]",1], res[u,:"K[6]",1] ]
    d = [ res[u,:"d[1]",1], res[u,:"d[2]",1], res[u,:"d[3]",1], res[u,:"d[4]",1], res[u,:"d[5]",1], res[u,:"d[6]",1] ]
    v = [ res[u,:"v[1]",1], res[u,:"v[2]",1], res[u,:"v[3]",1], res[u,:"v[4]",1] ]
    p = ( b, K, d, v, tau, hsa )
     
    msol = solve( 
      remake( prob, u0=u0, h=h, tspan=tspn, p=p  ), 
      solver, 
      callback=cb,  
      saveat=dt   ) #
    # plot!( msol; alpha=0.05 )
    plot!( msol.t, reduce(hcat, msol.u)'[:,j], color=[2 2], alpha=0.05 ) 

  end
end

plot!(; legend=false, xlim=(1997,2023) )


u0 = [ 
  mean( res[[:"K[1]"]].value ),
  mean( res[[:"K[2]"]].value ),
  mean( res[[:"K[3]"]].value ),
  mean( res[[:"K[4]"]].value ),
  mean( res[[:"K[5]"]].value ),
  mean( res[[:"K[6]"]].value )
] .*  [ 
  mean( res[[:"m1[1]"]].value ),
  mean( res[[:"m2[1]"]].value ),
  mean( res[[:"m3[1]"]].value ),
  mean( res[[:"m4[1]"]].value ),
  mean( res[[:"m5[1]"]].value ),
  mean( res[[:"m6[1]"]].value )
]



b = [ mean( res[[:"b[1]"]].value), mean( res[[:"b[2]"]].value) ]
K = [ mean( res[[:"K[1]"]].value), mean( res[[:"K[2]"]].value), 
      mean( res[[:"K[3]"]].value), mean( res[[:"K[4]"]].value),
      mean( res[[:"K[5]"]].value), mean( res[[:"K[6]"]].value) ]  ; 
d = [ mean( res[[:"d[1]"]].value), mean( res[[:"d[2]"]].value), 
      mean( res[[:"d[3]"]].value), mean( res[[:"d[4]"]].value),
      mean( res[[:"d[5]"]].value), mean( res[[:"d[6]"]].value) ]   
v = [ mean( res[[:"v[1]"]].value), mean( res[[:"v[2]"]].value), 
      mean( res[[:"v[3]"]].value), mean( res[[:"v[4]"]].value) ]

pm = ( b, K, d, v, tau, hsa ) 


msol = solve( remake( prob, u0=u0, h=h, tspan=tspan, p=pm ), solver, callback=cb, saveat=dt )  
plot!(msol, label="ode-mean-fishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 
plot!(; legend=false, xlim=(1997,2023) )

prob2 = DDEProblem( size_structured!, u0, h, tspan, pm, saveat=dt )
msol = solve( prob2,  solver, saveat=dt ) #  effective nullify callbacks
plot!(msol, label="ode-mean-nofishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 

# back transform S1 to normal scale 
yhat = ( S1 .* mean(res[[:"q[1]"]].value) .- mean(res[[:"qc[1]"]].value )) .* mean(res[[:"K[1]"]].value) 
scatter!(survey_time, yhat   ; color=[1 2])
plot!(survey_time, yhat  ; color=[1 2])
plot!(; legend=false, xlim=(1997,2023) )


# sample and plot means from model
j = 1  # state variable index
j = 6
w = zeros(nT)
for u in 1:length(res)  
  for i in 1:nT
    w[i] = res[u,Symbol("K[$j]"),1] * res[u,Symbol("m$j[$i]"),1]
  end
  plot!(survey_time, w  ;  alpha=0.1, color=[j j])
end

plot!(; legend=false, xlim=(1997,2023) )

u = zeros(nT)
v = zeros(nT)
for  i in 1:nT
  u[i] = mean( res[:,Symbol("m$j[$i]"),:] .* res[:,Symbol("K[$j]"),:] ) 
  v[i] = std(  res[:,Symbol("m$j[$i]"),:] .* res[:,Symbol("K[$j]"),:] ) 
end
scatter!(survey_time, u  ; color=[3 2])

  
# misc computed quantities

# params need to be named  .. return only that which is specified by "return()", below
pm = ( b=b, K=K, d=d, v=v, tau=tau, hsa=hsa ) 

@model function fm_test( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, N=length(S1), ::Type{T}=Float64 ) where {T}  
  
  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(N+M)
  B = zeros(N+M)
  C = zeros(N+M)

  C[1:N] = removed ./ K
  C[(N+1):(M+N)] = er .* bm[(N):(M+N-1)]
  C = 1.0 .- C / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = bm .* K  
 
  # recaled estimates
  B[1:N] = bm[1:N] *. K - L[1:N] ;
  B[(N+1):(M+N)] = (bm[(N+1):(M+N)] - C[(N):(M+N-1)]) *. K ;
 
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




# ------------------------------
# Delay SSA size structured

# https://github.com/palmtree2013/DelaySSAToolkit.jl
# https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/
# NOTE::: require 03.snowcrab_carstm.r to be completed 


# from https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/:

# "A major assumption behind the majority of stochastic models of biochemical kinetics is the memoryless hypothesis, i.e., the stochastic dynamics of the reactants is only influenced by the current state of the system, which implies that the waiting times for reaction events obey exponential distributions. Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics for such systems [1]. While this Markovian assumption considerably simplifies model analysis, it is dubious for modelling certain non-elementary reaction events that encapsulate multiple intermediate reaction steps [2]."

# For a few number of jumps, DelayRejection and DelayDirect will often perform better than other aggregators.

# For large numbers of jumps with sparse chain like structures and similar jump rates, for example continuous time random walks, DelayDirectCR and DelayMNRM often have the best performance.

 
# NOTE::: require 03.snowcrab_carstm.r to be completed 


dir = expanduser("~/julia/snowcrab/")  # The directory of your package, for you maybe "C:\something"  
push!(LOAD_PATH, dir)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(dir)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions", "DynamicPPL",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  "Interpolations",
  "Plots", "StatsPlots", "MultivariateStats", "Graphviz_jll", "DelaySSAToolkit", "JumpProcesses"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages



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
     
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
    o = load( fndat, convert=true)
    Y = o["Y"]
    Kmu = o["Kmu"]
    removals = o["L"]

end

plot(  Y[:,:yr], Y[:,:cfasouth_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfasouth_M1] )
plot!( Y[:,:yr] .+2, Y[:,:cfasouth_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfasouth_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfasouth_M4] )


plot(  Y[:,:yr], Y[:,:cfanorth_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfanorth_M1] )
plot!( Y[:,:yr] .+2 , Y[:,:cfanorth_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfanorth_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfanorth_M4] )


plot(  Y[:,:yr], Y[:,:cfa4x_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfa4x_M1] )
plot!( Y[:,:yr] .+2 , Y[:,:cfa4x_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfa4x_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfa4x_M4] )



# ------------------------------

Turing.setprogress!(false);
# Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
# Turing.setadbackend(:tracker)
 
rn = @reaction_network begin
  v1 * S2,             S2 --> S1
  v2 * S3,             S3 --> S2
  v3 * S4,             S4 --> S3
  v4 * S5,             S5 --> S4
  b1 * S6,             S6 --> S5
  b2 * S6,             S6 --> S6
end v1 v2 v3 v4 b1 b2 

species(rn)
parameters(rn)
Graph(rn)


tspan = (0.0, 100.0)
dt = 0.1

u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 100

b=[ 3.0, 1.24]
K=[100.0, 100.0, 100.0, 100.0, 100.0, 100.0];
d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
v=[0.8, 0.8, 0.8, 0.8];  
p = ( 
  b1=b[1], b2=b[2],
  K1=K[1], K2=K[2], K3=K[3], K4=K[4], K5=K[5], K6=K[6], 
  d1=d[1], d2=d[2], d3=d[3], d4=d[4], d5=d[5], d6=d[6], 
  v1=v[1], v2=v[2], v3=v[3], v4=v[4] 
)

# instantaneous probability per time a jump occurs when
# the current state is u, current parameters are p, and the time is t
b1_rate(u,p,t) = p.b1
b1_affect!(integrator) = (integrator.u[1] += 1)
b1_jump = ConstantRateJump(b1_rate, b1_affect!)

d1_rate(u,p,t) = p.d1 * u[1]
d1_affect!(integrator) = (integrator.u[1] -= 1; integrator.u[2] += 1)
d1_jump = ConstantRateJump(d1_rate, d1_affect!)

h1_rate(u,p,t) = hsa(t,1) 
h1_affect!(integrator) = (integrator.u[1] -= 1);
h1_jump = VariableRateJump(h1_rate, h1_affect!)

# this is fora comre complex random rate (nonconstant rate) jump .. 
rng = JumpProcesses.DEFAULT_RNG
Random.TaskLocalRNG()
b1r_rate(u,p,t) = p.b1
# define the affect function via a closure

b1r_affect! = integrator -> let rng=rng
    # N(t) <-- N(t) + 1
    integrator.u[1] += 1
    # G(t) <-- G(t) + C_{N(t)}
    integrator.u[2] += rand(rng, (-1,1))
    nothing
end
b1r_jump = ConstantRateJump(b1r_rate, b1r_affect!)


dprob = DiscreteProblem( u0, tspan, p )
jprob = JumpProblem(dprob, Direct(), b1_jump, h1_jump, d1_jump; save_positions = (false, false) )
sol = solve(jprob, SSAStepper(); saveat=dt )
plot(sol, label="N(t)", xlabel="t", legend=:bottomright)

 

delay_trigger_affect1! = function (integrator, rng)
  append!(integrator.de_chan[1], 1)
end
delay_trigger_affect2! = function (integrator, rng)
  append!(integrator.de_chan[2], 1)
end
delay_trigger_affect3! = function (integrator, rng)
  append!(integrator.de_chan[3], 1)
end
delay_trigger_affect4! = function (integrator, rng)
  append!(integrator.de_chan[4], 1)
end
delay_trigger_affect5! = function (integrator, rng)
  append!(integrator.de_chan[5], 8)
end
delay_trigger_affect6! = function (integrator, rng)
  append!(integrator.de_chan[6], 8)
end
 
delay_trigger = Dict(
  7=>delay_trigger_affect1!,
  8=>delay_trigger_affect2!,
  9=>delay_trigger_affect3!,
  10=>delay_trigger_affect4!,
  11=>delay_trigger_affect5!,
  12=>delay_trigger_affect6!
)

delay_interrupt = Dict()
delay_complete = Dict(
  1=>[1=>1, 2=>-1],
  2=>[2=>1, 3=>-1],
  3=>[3=>1, 4=>-1],
  4=>[4=>1, 5=>-1],
  5=>[5=>1],
  6=>[6=>1]
)

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

de_chan0 = [[]]

djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0 )
  
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, 
  saveat=dt, save_positions=(false,false))

sol = solve(djprob, SSAStepper())


function size_structured!( du, u, h, p, t)
  # S1, S2, S3, S4, S5, S6 = u
  b, K, d, v, tau, hsa  = p
  tr21 = v[1] * h(p, t-1)[2]   # transition 2 -> 1   
  tr32 = v[2] * h(p, t-1)[3]   # transitiom 3 -> 2
  tr43 = v[3] * h(p, t-1)[4]   # transitiom 4 -> 3
  tr54 = v[4] * h(p, t-1)[5]   # transitiom 5 -> 4
  f8  = h(p, t-8)[6]           # no fem 8 yrs ago
  du[1] = tr21             - (d[1] * u[1]) * (u[1]/ K[1]) * hsa(t,1)       
  du[2] = tr32      - tr21 - (d[2] * u[2]) * (u[2]/ K[2]) * hsa(t,2) 
  du[3] = tr43      - tr32 - (d[3] * u[3]) * (u[3]/ K[3]) * hsa(t,3)
  du[4] = tr54      - tr43 - (d[4] * u[4]) * (u[4]/ K[4]) * hsa(t,4)
  du[5] = b[1] * f8 - tr54 - (d[5] * u[5]) * (u[5]/ K[5]) * hsa(t,5) 
  du[6] = b[2] * f8        - (d[6] * u[6]) * (u[6]/ K[6]) * hsa(t,6)  # fem mat simple logistic with lag tau and density dep on present numbers
end




if false
  # testing DDE version -- initial values for size_structured! 
  
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 100
  b=[ 3.0, 1.24]
  K=[100.0, 100.0, 100.0, 100.0, 100.0, 100.0];
  d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
  v=[0.8, 0.8, 0.8, 0.8];  
  tau=1.0 

  ps = [ d, K, v, b ]

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
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
  tau = 1  # delay
  lags = [tau]
  solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
  prob = DDEProblem( size_structured , u0, h, tspan, p; constant_lags=lags )
  res =  solve( prob,  solver, saveat=0.1 )
  plot( res )


  # ---------------
  #test runs

  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 1.e9
  b=[3.0, 1.24]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*1.e-9; 
  d=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
  v=[0.9, 0.9, 0.9, 0.8];  
  tau=1.0; 

  p = ( b, K, d, v, tau, hsa )   

  forcing_time = survey_time

  external_forcing = ones(length(forcing_time),6)  # turns it off
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  plot!( msol, label="dde, no hsa, no fishing", legend=:left )
   

  external_forcing = rand(length(forcing_time),6) # random
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  plot!( msol, label="dde, with hsa, no fishing", legend=:left )


end




# -------------------------
# other parameters

au = 2  # cfasouth
aulab ="cfasouth"
eps = 1.0e-9

# convert to number .. 0.56 is ave mean weight of fb
kmu = Kmu[au] * 1000 *1000 / 0.56

# "survey index"
S1 = Y[:,Symbol("$aulab","_M0")]  
S2 = Y[:,Symbol("$aulab","_M1")]  
S3 = Y[:,Symbol("$aulab","_M2")]  
S4 = Y[:,Symbol("$aulab","_M3")] 
S5 = Y[:,Symbol("$aulab","_M4")]  
S6 = Y[:,Symbol("$aulab","_f_mat")]  

nS = 6 # n components
nT = length(S1)
dt = 0.1
yrs = 1999:2021
tspan = (1990.0, 2030.0)

survey_time = Y[:,:yrs]   # time of observations for survey

# this only adds habitat space  ... predation is also a useful one .. 
# speed is the issue 
forcing_time = survey_time

# invert sa to fraction

external_forcing = reshape( [
    1.0 .- Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
    1.0 .- Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] ) 
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
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5 *kmu
 
tau = 1  # delay
lags = [tau]

solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()

 
 
# ---------------


@model function fishery_model_turing_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, nS=6, N=length(S1), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K ~ filldist( truncated( Cauchy( kmu, kmu*0.25), kmu/10.0, kmu*10.0), nS )  

    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ filldist( truncated( Cauchy( 1.0, 1.0), 0.1, 10.0), nS )  
    qc ~ filldist( truncated( Cauchy(0.0, 0.1, -1.0, 1.0), nS )  
  
    # initial conditions
    m1 =  Vector{T}(undef, N)
    m2 =  Vector{T}(undef, N)
    m3 =  Vector{T}(undef, N)
    m4 =  Vector{T}(undef, N)
    m5 =  Vector{T}(undef, N)
    m6 =  Vector{T}(undef, N)

    m1[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m2[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m3[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m4[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m5[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event
    m6[1] ~  truncated( Cauchy( 0.5, 0.5), 0.1, 1.0 )  ; # starting b prior to first catch event

    # birth rate from F_8
    b ~ filldist( truncated( Cauchy(1.0, 0.1), 0.25, 3.0), 2 ) 
     
    # mortality
    d ~ filldist( truncated( Cauchy(0.2, 0.1), 1.0e-9, 0.9), nS )  

    # transition rates
    v ~ filldist( truncated( Cauchy(0.9, 0.1), 1.0e-9, 1.0 ), 4 ) 
 
    # process model
    u0 = T[ m1[1], m2[1], m3[1], m4[1], m5[1], m6[1] ] .* K 
    p = ( b, K, d, v, tau, hsa )
    msol = solve( remake( prob; u0=u0, h=h, tspan=tspan, p=p ), solver, callback=cb, saveat=dt )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        usk = max.( msol.u[j[1]] ./ K, 1.0e-9 )
        m1[i] ~ TruncatedNormal( usk[1], bpsd, 1.0e-9, 1.0)  ; 
        m2[i] ~ TruncatedNormal( usk[2], bpsd, 1.0e-9, 1.0)  ; 
        m3[i] ~ TruncatedNormal( usk[3], bpsd, 1.0e-9, 1.0)  ; 
        m4[i] ~ TruncatedNormal( usk[4], bpsd, 1.0e-9, 1.0)  ; 
        m5[i] ~ TruncatedNormal( usk[5], bpsd, 1.0e-9, 1.0)  ; 
        m6[i] ~ TruncatedNormal( usk[6], bpsd, 1.0e-9, 1.0)  ; 
      end
    end
  
    # observation model
    @. S1 ~ TruncatedNormal( (m1 + qc[1]) * q[1], bosd, 0.0, 1.0 )   
    @. S2 ~ TruncatedNormal( (m2 + qc[2]) * q[2], bosd, 0.0, 1.0 )   
    @. S3 ~ TruncatedNormal( (m3 + qc[3]) * q[3], bosd, 0.0, 1.0 )   
    @. S4 ~ TruncatedNormal( (m4 + qc[4]) * q[4], bosd, 0.0, 1.0 )   
    @. S5 ~ TruncatedNormal( (m5 + qc[5]) * q[5], bosd, 0.0, 1.0 )   
    @. S6 ~ TruncatedNormal( (m6 + qc[6]) * q[6], bosd, 0.0, 1.0 )   
    
end


# ---------------

@model function fishery_model_turing_incremental_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, N=length(S1), ::Type{T}=Float64 ) where {T}  
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors
  K ~ filldist( truncated( Cauchy( kmu, kmu*0.25), kmu/10.0, kmu*10.0), nS )  

  bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ filldist( truncated( Cauchy( 1.0, 1.0), 0.1, 10.0), nS )  
  qc ~ filldist( truncated( Cauchy(0.0, 0.1, -1.0, 1.0), nS )  
 
   # initial conditions
  m1 =  Vector{T}(undef, N)
  m2 =  Vector{T}(undef, N)
  m3 =  Vector{T}(undef, N)
  m4 =  Vector{T}(undef, N)
  m5 =  Vector{T}(undef, N)
  m6 =  Vector{T}(undef, N)

  m1[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
  m2[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
  m3[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
  m4[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
  m5[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
  m6[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event

 
  # birth rate from F_8
  b ~ filldist( truncated( Cauchy(1.0, 0.1), 0.25, 3.0), 2 ) 
   
  # mortality
  d ~ filldist( truncated( Cauchy(0.2, 0.1), 1.0e-9, 0.9), nS )  

  # transition rates
  v ~ filldist( truncated( Cauchy(0.9, 0.1), 1.0e-9, 1.0 ), 4 ) 
 

  # process model
  t0 = floor(survey_time[1]) 
  
  for i in 2:N
    tsp = (t0+i-1.1, t0+i+0.1 )
    u0 = T[ m1[1], m2[1], m3[1], m4[1], m5[1], m6[1] ] .* K 
    p = ( b, K, d, v, tau, hsa )
  
    msol = solve( 
      remake( prob; u0=u0, h=h, tspan=tsp, p=p, constant_lags=lags   ), 
      solver, 
      callback=cb,  
      saveat=dt 
    )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    j = findall(t -> t==survey_time[i], msol.t)
    if length(j) > 0
      usk = max.( msol.u[j[1]] ./ K, 1.0e-9 )
      m1[i] ~ TruncatedNormal( usk[1], bpsd, 1.0e-9, 1.0)  ; 
      m2[i] ~ TruncatedNormal( usk[2], bpsd, 1.0e-9, 1.0)  ; 
      m3[i] ~ TruncatedNormal( usk[3], bpsd, 1.0e-9, 1.0)  ; 
      m4[i] ~ TruncatedNormal( usk[4], bpsd, 1.0e-9, 1.0)  ; 
      m5[i] ~ TruncatedNormal( usk[5], bpsd, 1.0e-9, 1.0)  ; 
      m6[i] ~ TruncatedNormal( usk[6], bpsd, 1.0e-9, 1.0)  ; 
    end
  end

  # observation model
  @. S1 ~ TruncatedNormal( (m1 + qc[1]) * q[1], bosd, 0.0, 1.0 )   
  @. S2 ~ TruncatedNormal( (m2 + qc[2]) * q[2], bosd, 0.0, 1.0 )   
  @. S3 ~ TruncatedNormal( (m3 + qc[3]) * q[3], bosd, 0.0, 1.0 )   
  @. S4 ~ TruncatedNormal( (m4 + qc[4]) * q[4], bosd, 0.0, 1.0 )   
  @. S5 ~ TruncatedNormal( (m5 + qc[5]) * q[5], bosd, 0.0, 1.0 )   
  @. S6 ~ TruncatedNormal( (m6 + qc[6]) * q[6], bosd, 0.0, 1.0 )   
  
end

# ---------------

  
prob = DDEProblem( size_structured!, u0, h, tspan, p, constant_lags=lags )
fmod = fishery_model_turing_dde( S1, S2, S3, S4, S5, S6, kmu, tspan, prob )
fmod = fishery_model_turing_incremental_dde( S1, S2, S3, S4, S5, S6, kmu,  tspan, prob )
 

# testing
res  =  sample( fmod,  Turing.MH(), 3 )
res  =  sample( fmod,  Turing.NUTS( 3, 0.65), 3 )


n_samples = 1000
n_adapts = 500

res  =  sample( fmod,  Turing.NUTS(n_adapts, 0.65), n_samples )

t0 = floor(survey_time[1])

v = 1 # S1
v = 6 # female

for u in 1:100
  for  i in 1:N
    u0 = [ 
      res[u,Symbol("m1[$i]"),1] * res[u,:K1,1],
      res[u,Symbol("m2[$i]"),1] * res[u,:K2,1],
      res[u,Symbol("m3[$i]"),1] * res[u,:K3,1],
      res[u,Symbol("m4[$i]"),1] * res[u,:K4,1],
      res[u,Symbol("m5[$i]"),1] * res[u,:K5,1],
      res[u,Symbol("m6[$i]"),1] * res[u,:K6,1]
    ]
    
    tspn = ( t0+i-1.1, t0+i+0.1 )

    p = ( 
      res[u,:b4,1], res[u,:b6,1], 
      res[u,:K1,1], res[u,:K2,1], res[u,:K3,1], res[u,:K4,1], res[u,:K5,1], res[u,:K6,1], 
      res[u,:d1,1], res[u,:d2,1], res[u,:d3,1], res[u,:d4,1], res[u,:d5,1], res[u,:d6,1],
      res[u,:v1,1], res[u,:v2,1], res[u,:v3,1], res[u,:v4,1], 
      tau, hsa
    )
  
    msol = solve( 
      remake( prob, u0=u0, h=h, tspan=tspn, p=p  ), 
      solver, 
      callback=cb,  
      saveat=dt   ) #
    # plot!( msol; alpha=0.05, color=[3 4] )
    plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[2 2], alpha=0.05 ) 

  end
end


plot!(; legend=false, xlim=(1997,2023) )

u0 = [ 
  mean( res[[:"K[1]"]].value ),
  mean( res[[:"K[2]"]].value ),
  mean( res[[:"K[3]"]].value ),
  mean( res[[:"K[4]"]].value ),
  mean( res[[:"K[5]"]].value ),
  mean( res[[:"K[6]"]].value )
] .*  [ 
  mean( res[[:"m1[1]"]].value ),
  mean( res[[:"m2[1]"]].value ),
  mean( res[[:"m3[1]"]].value ),
  mean( res[[:"m4[1]"]].value ),
  mean( res[[:"m5[1]"]].value ),
  mean( res[[:"m6[1]"]].value )
]



b = [ mean( res[[:"b[1]"]].value), mean( res[[:"b[2]"]].value) ]
K = [ mean( res[[:"K[1]"]].value), mean( res[[:"K[2]"]].value), 
      mean( res[[:"K[3]"]].value), mean( res[[:"K[4]"]].value),
      mean( res[[:"K[5]"]].value), mean( res[[:"K[6]"]].value) ]  ; 
d = [ mean( res[[:"d[1]"]].value), mean( res[[:"d[2]"]].value), 
      mean( res[[:"d[3]"]].value), mean( res[[:"d[4]"]].value),
      mean( res[[:"d[5]"]].value), mean( res[[:"d[6]"]].value) ]   
v = [ mean( res[[:"v[1]"]].value), mean( res[[:"v[2]"]].value), 
      mean( res[[:"v[3]"]].value), mean( res[[:"v[4]"]].value) ]

pm = ( b, K, d, v, tau, hsa ) 


msol = solve( remake( prob, u0=u0, h=h, tspan=tspan, p=pm ), solver, callback=cb, saveat=dt )  
plot!(msol, label="ode-mean-fishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 
plot!(; legend=false, xlim=(1997,2023) )

prob2 = DDEProblem( size_structured!, u0, h, tspan, pm, saveat=dt )
msol = solve( prob2,  solver, saveat=dt ) #  effective nullify callbacks
plot!(msol, label="ode-mean-nofishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 

# back transform S1 to normal scale 
yhat = ( S1 .* mean(res[[:q1]].value) .- mean(res[[:qc1]].value )) .* mean(res[[:K1]].value) 
scatter!(survey_time, yhat   ; color=[1 2])
plot!(survey_time, yhat  ; color=[1 2])
plot!(; legend=false, xlim=(1997,2023) )


# sample and plot means from model
j = 1  # state variable index

w = zeros(N)
for u in 1:length(res)  
  for i in 1:N
    w[i] = res[u,Symbol("K[$j]"),1] * res[u,Symbol("m$j[$i]"),1]
  end
  plot!(survey_time, w  ;  alpha=0.1, color=[j j])
end

plot!(; legend=false, xlim=(1997,2023) )

u = zeros(N)
v = zeros(N)
for  i in 1:N
  u[i] = mean( res[:,Symbol("m$j[$i]"),:] .* res[:,Symbol("K[$j]"),:] ) 
  v[i] = std(  res[:,Symbol("m$j[$i]"),:] .* res[:,Symbol("K[$j]"),:] ) 
end
scatter!(survey_time, u  ; color=[3 2])


# misc computed quantities

# params need to be named  .. return only that which is specified by "return()", below
pm = ( b=b, K=K, d=d, v=v, tau=tau, hsa=hsa ) 

@model function fm_test( S1, S2, S3, S4, S5, S6, kmu, tspan, prob, N=length(S1), ::Type{T}=Float64 ) where {T}  
  # deterministic computations: do from similations:
  M = 5
  er = 0.2
  eps = 1.0e-9
  
  F = zeros(N+M)
  B = zeros(N+M)
  C = zeros(N+M)

  C[1:N] = removed ./ K
  C[(N+1):(M+N)] = er .* m1[ (N):(M+N-1)]
  C = 1.0 .- C / m1

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r * exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = m1 .* K  
 
  # recaled estimates
  B[1:N] = m1[1:N] *. K - L[1:N] ;
  B[(N+1):(M+N)] = (m1[(N+1):(M+N)] - C[(N):(M+N-1)]) *. K ;
 
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


