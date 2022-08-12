
# ----------------------------------------------
# ODE

# project_directory = string(expanduser("~/projects/dynamical_model/"), "snowcrab")
# # project_directory = @__DIR__() #  same folder as the file

# push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

# import Pkg  # or using Pkg
# Pkg.activate(project_directory)  # so now you activate the package
# # Pkg.activate(@__DIR__()) #  same folder as the file itself.

# Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  
  "Plots", "StatsPlots", "MultivariateStats"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages



# ----------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

# NOTE::: require 03.snowcrab_carstm.r to be completed 
 
 
get_data_with_RCall = false

if get_data_with_RCall

    using RCall

    # typing <$> in Julia's  command prompt starts an R session.  
    
    $
      
    {
        # this is R-code that creates local RData file with required data
        source( file.path( code_root, "bio_startup.R" )  )
        require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        fishery_model_data_inputs( year.assessment=2021, type="numerical_dynamics" )
        fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics" )
        fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics" ) 
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


# ----------------------------------------------


Turing.setprogress!(false);
Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
# Turing.setadbackend(:tracker)
 
 
 
function Logistic!( du, u, p, t)
  du[1] = p[1] * u[1]  * (1.0 - u[1]/p[2])  # specific rate
end


if false
  # testing ODE version -- initial values for Logistic! 
  u0 = [0.1 ]
  p = (1.0, 1.0)
  tspan = (0.0, 5.0)
  prob = ODEProblem( Logistic!, u0,  tspan, p )
  res =  solve( prob, Tsit5(), saveat=0.1 )
  plot!( res )
  # plot( res.t, reduce(hcat, res.u)' )  #res.u is an Vector{Vector{}} 
end

 

# -------------------------
# other parameters

au = 2  # cfasouth
eps = 1e-9

tspan = (1999.0, 2025.0)

# convert to number .. 0.56 is ave mean weight of fb
kmu = Kmu[au] * 1000 *1000 / 0.56
ksd = kmu * 0.2
  
M0 = Y[:,:cfasouth_M0]  # "survey index"
survey_time = Y[:,:yrs]  # time of observations for survey

N = length(M0)
dt = 0.1

fish_time = removals[:,:ts]
removed = removals[:,:cfasouth]

function affect_fishing!(integrator)
  i = findall(t -> t==integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
end

cb =  PresetTimeCallback( fish_time, affect_fishing! )
 
p = ( 1.0, kmu, rand(Beta(2,1))*kmu, tspan[1], tspan[2] )  #p[4] tspan[0], p[5] is dtspan



# ---------------
if false
    #test run
    
    prob = ODEProblem( Logistic!, [0.1], tspan, p )
    msol =  solve( prob, Tsit5(), saveat=dt )
    plot( msol, label="ode, no fishing", legend=:left )
    
    # test 2 .. alt y0
    prob = ODEProblem( Logistic!, [p[3]], tspan_func, p )
    msol =  solve( prob, Tsit5(), saveat=dt )
    plot!( msol, label="ode, no fishing random start 1" )
    
    # test 3
    prob = remake(prob; f=Logistic!, u0=[0.1], tspan=tspan_func(p), p=p )
    msol =  solve( prob, Tsit5(), saveat=dt )
    plot!( msol, label="ode, no fishing random start 2" )
    
    #test 4 with removals
    prob = remake(prob; f=Logistic!, u0=[0.1], tspan=tspan_func(p), p=p, callback=cb )
    msol =  solve( prob, Tsit5(), saveat=dt  )
    plot!( msol, label="ode, fishing", legend = :bottomright, ylim=(0,kmu)) 
    
end



# ---------------

@model function fishery_model_turing_ode( M0, kmu, ksd, tspan, prob, N=length(M0), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K  ~  TruncatedNormal( kmu, ksd, kmu/5.0, kmu*5.0)   ; # (mu, sd)
    r ~  TruncatedNormal( 1.0, 0.25, 0.25, 2.0)   # (mu, sd)
    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    q ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # initial conditions
    m0 =  Vector{T}(undef, N)
    m0[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event

    # process model
    u0 = T[m0[1]*K]
    p = [ r, K ]
    msol = solve( 
      remake( prob, f=Logistic!, u0=u0, tspan=tspan, p=p ), 
      Rosenbrock23(), 
      callback=cb,  
      saveat=dt 
    ) 
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        ym = msol.u[j[1]][1] / K
        if typeof( ym ) !== Float64
         # Turing.@addlogprob! -Inf
         #  return nothing
           ym = m0[i-1] 
        end
        m0[i] ~ TruncatedNormal( max( ym, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
      end
    end
 
    # observation model
    # @. M0 ~ TruncatedNormal( (m0 *q) + qc, bosd, -1.0, 1.0 )  # M0 in SD units 
    @. M0 ~ TruncatedNormal( (m0 + qc) * q, bosd, -1.0, 1.0 )  # M0 in SD units 

end



@model function fishery_model_turing_incremental_ode( M0, kmu, ksd, tspan, prob, N=length(M0), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K  ~  TruncatedNormal( kmu, ksd, kmu/5.0, kmu*5.0)   ; # (mu, sd)
    r ~  TruncatedNormal( 1.0, 0.25, 0.25, 2.0)   # (mu, sd)
    bpsd ~  truncated( Cauchy( 0.0, 1.0), 1e-9, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  truncated( Cauchy( 0.0, 1.0), 1e-9, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    q  ~  truncated( Cauchy( 0.0, 1.0), 1e-9, 5.0 )   ; # i.e., Y:b scaling coeeficient
    qc ~  truncated( Cauchy( 0.0, 1.0), -1.0, 1.0 ) ; # i.e., Y:b offset constant   
    
    # initial conditions
    m0 =  Vector{T}(undef, N)
    # m0 =tzeros(N)
    m0[1] ~ truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    t0 = floor(survey_time[1])
    
    # process model
    for i in 2:N
      u0 = T[m0[i-1]*K]
      tsp = (t0+i-1.1, t0+i+0.1 )
      p = [r, K]
      prob2 = remake(prob; f=Logistic!, u0=u0, tspan=tsp, p=p )
      msol = solve( prob2, Rosenbrock23(), callback=cb, saveat=dt ) 
      if msol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
      end
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        m0[i] ~ TruncatedNormal(   msol.u[j[1]][1] / K, bpsd, 1e-9, 1.25)  ; 
      end
    end

    # observation model
    # @. M0 ~ Cauchy( (m0 .+ qc) .* q, bosd  ) 
    @. M0 ~ TruncatedNormal( (m0 + qc) * q, bosd, -1.0, 1.0 )  # M0 in SD units 

end

  

#  ----
#  run
 
prob = ODEProblem( Logistic!, [rand(Beta(5,1))*kmu], tspan, [1.0, kmu], saveat=dt, callback=cb    )
 
fmod = fishery_model_turing_incremental_ode( M0, kmu, ksd,  tspan, prob )

fmod = fishery_model_turing_ode( M0, kmu, ksd,  tspan, prob )

 
# testing
n_samples = 3
n_adapts = 3
n_chains = 1


# production
n_samples = 1000
n_adapts = 500
n_chains = 4

sampler = Turing.MH()
sampler = Turing.NUTS(n_adapts, 0.65)

res  =  sample( fmod, sampler, MCMCThreads(), n_samples, n_chains )
# if on windows o threads not working:
# res = mapreduce(c -> sample(fmod, sampler, n_samples), chainscat, 1:n_chains)

show(stdout, "text/plain", summarize(res))
 
histogram(res[:p])


t0 = floor(survey_time[1])

for u in 1:1000 
  for  i in 1:N
    u0 = [res[u,:K,1] * res[u,Symbol("m0[$i]"),1]]
    tsp = ( t0+i-1.1, t0+i+0.1 )
  
    msol = solve( 
      remake( prob, f=Logistic!, u0=u0, tspan=tsp, p=[res[u,:r,1], res[u,:K,1]] ), 
      Rosenbrock23(), 
      callback=cb,  
      saveat=dt   ) #
    plot!(msol; alpha=0.05, color=[3 4], xlim=tsp)
  end
end

  plot!(; legend=false, ylim=(0,kmu*1.75) )

  k0 = [(floor(mean( res[[:"m0[1]"]].value ) *  mean( res[[:K]].value ) ))]

  pm = [mean( res[[:r]].value ), mean( res[[:K]].value ) ]

  msol = solve( remake( prob, u0=k0, tspan=tspan, p=pm ), Rosenbrock23(), callback=cb, saveat=dt )  
  plot!(msol, label="ode-mean-fishing")

  prob2 = ODEProblem( Logistic!, k0, tspan, pm, saveat=dt )
  msol = solve( prob2, Tsit5(), saveat=dt ) #  effective nullify callbacks
  plot!(msol, label="ode-mean-nofishing")

  # back transform M0 to normal scale 
  yhat = ( M0  .- mean(res[[:"qc"]].value)) ./ mean(res[[:"q"]].value) .* mean(res[[:"K"]].value) 
  scatter!(survey_time, yhat   ; color=[1 2])
  plot!(survey_time, yhat  ; color=[1 2])


  
  w = zeros(N)
  for u in 1:1000  
    for i in 1:N
      w[i] = res[u,:K,1] * res[u,Symbol("m0[$i]"),1]
    end
    plot!(survey_time, w  ;  alpha=0.1, color=[5 2])
  end

  u = zeros(N)
  v = zeros(N)

  for  i in 1:N
    u[i] = mean( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
    v[i] = std( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
  end
  scatter!(survey_time, u  ; color=[5 2])
  
  
  # look at predictions:
  si_pred = Vector{Union{Missing, Float64}}(undef, length(M0))
  fmod_pred = fmod( si_pred, kmu, ksd,  tspan, prob  ) 
 
  predictions = predict(fmod_pred, res)
  y_pred = vec(mean(Array(group(predictions, :M0)); dims = 1));
  
  plot( M0, y_pred )
  sum(abs2, M0 - y_pred) ≤ 0.1


 
 
 
 
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


   


  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(sN+M)
  B = zeros(N+M)
  C = zeros(N+M)

  C[1:N] = removals ./ K
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


 
using Flux, DiffEqFlux
params = Flux.params(p)

msol =  solve( prob, Tsit5(), callback=cb, saveat=dt)  

function predict_rd() # Our 1-layer "neural network"
  solve(prob,Tsit5(),p=p,saveat=dt)[1] # override with new parameters
end

loss_rd() = sum(abs2,x-1 for x in predict_rd()) # loss function (squared absolute) x-1

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cbflux = function () #callback
  # function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=p),Tsit5(),saveat=dt), ylim=(0,kmu*2)))
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
  x -> solve(prob,Tsit5(),u0=x,saveat=0.1)[1,:],
  Dense(288, 10), softmax) |> gpu


