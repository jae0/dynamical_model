# ----------------------------------------------
## SDE
# https://diffeq.sciml.ai/stable/solvers/sde_solve/
# Recommended Methods

# For most Ito diagonal and scalar noise problems where a good amount of accuracy is required and mild stiffness may be an issue, the SOSRI algorithm should do well. If the problem has additive noise, then SOSRA will be the optimal algorithm. At low tolerances (<1e-4?) SRA3 will be more efficient, though SOSRA is more robust to stiffness. For commutative noise, RKMilCommute is a strong order 1.0 method which utilizes the commutivity property to greatly speed up the stochastic iterated integral approximation and can choose between Ito and Stratonovich. For non-commutative noise, difficult problems usually require adaptive time stepping in order to be efficient. In this case, LambaEM and LambaEulerHeun are adaptive and handle general non-diagonal problems (for Ito and Stratonovich interpretations respectively). If adaptivity isn't necessary, the EM and EulerHeun are good choices (for Ito and Stratonovich interpretations respectively).

# # For stiff problems with additive noise, the high order adaptive method SKenCarp is highly preferred and will solve problems with similar efficiency as ODEs. If possible, stiff problems should be converted to make use of this additive noise solver. If the noise term is large/stiff, then the split-step methods are required in order for the implicit methods to be stable. For Ito in this case, use ISSEM and for Stratonovich use ISSEulerHeun. These two methods can handle any noise form.
# project_directory = string(expanduser("~/projects/dynamical_model/"), "snowcrab")
# # project_directory = @__DIR__() #  same folder as the file

# push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found

# import Pkg  # or using Pkg
# Pkg.activate(project_directory)  # so now you activate the package
# # Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",  "Arpack",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions", 
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  
  "Plots", "StatsPlots", "MultivariateStats"
]

#  Pkg.add( pkgs ) # add required packages

for pk in pkgs; @eval using $(Symbol(pk)); end
 

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
    include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_sde.jl")

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

# Turing.setprogress!(false);
# Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
 
 
rs = @reaction_network begin
  r, X--> 2X
  r*X/K, X --> 0
end r K 



if false

    # testing 
    u0 = [1]
    p = [1.0, 100]
    tspan = (0.0,10.0)
    dt = 0.1
        
    # basic dynamics  
    # (Direct(), DirectFW(), DirectCR(), SortingDirect(), RSSA(), FRM(), FRMFW(), NRM(), RSSACR(), RDirect())
    jprob = JumpProblem(rs, DiscreteProblem(rs, u0, tspan, p ), Direct())
    jsol = solve( jprob, SSAStepper(), saveat=dt  )
    plot!(jsol, lw=2, title=("Gillespie solution for $p"), label="simple"  )
    
    # as ode
    oprob = ODEProblem(rs, u0, tspan, p)
    osol =  solve( oprob, Tsit5(), saveat=0.1 )
    plot!( osol )

    # as an SPDE:
    @show rs_spde = convert(SDESystem, rs)
    prob = SDEProblem(rs_spde, u0, tspan, p)
    msol  = solve(prob, LambaEM(), reltol=1e-3)
    plot!(msol, lw=2, title="Adaptive SDE")
    
    
    Graph(rs)
    species(rs)
    parameters(rs)
    reactions(rs)
    latexify(rs, starred=true)
    latexify( )
    g = Graph(rs)
    
    using Plots, GraphRecipes
    plot(TreePlot(rs), method=:tree, fontsize=12, nodeshape=:ellipse)
 
end

 

# -------------------------
# other parameters

au = 2
eps = 1e-9
    
 
tspan = (1999.0, 2025.0)

# convert to number .. 0.56 is ave mean weight of fb
scalefactor = 1000
kmu = Kmu[au] * 1000 *1000 / 0.56  / scalefactor
ksd = Ksd[au] * 1000 *1000 / 0.56  / scalefactor

M0 = Y[:,:cfasouth]  # "survey index"
survey_time = Y[:,:yrs]  # time of observations for survey

N = length(M0)
dt = 0.1 

fish_time = removals[:,:ts]
removed = removals[:,:cfasouth] / scalefactor

function affect_fishing!(integrator)
  i = findall(t -> t==integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
end

cb =  PresetTimeCallback( fish_time, affect_fishing! )
 
 
p = ( 1.0, kmu )
u0 = [rand(Beta(2,1))*kmu]


 
# ---------------
if false
    #test run
  
  
    # simple run
    @show rs_spde = convert(SDESystem, rs)
    
    prob = SDEProblem(rs_spde, u0, tspan, p)
    msol  = solve(prob, LambaEM(), reltol=1e-3)
    plot!(msol, lw=2, title="Adaptive SDE")
      
    # simple run with fishing
    prob = SDEProblem(rs_spde, u0, tspan, p)
    msol  = solve(prob, LambaEM(), reltol=1e-3,  saveat=dt, callback=cb)
    plot!(msol, lw=2, title="Adaptive SDE", ylim=(0, kmu*1.5) )
     
    # testing use of remake
    u0 = [rand(Beta(5,1))*kmu]
    p = ( 1.0, kmu  )  #p[4] tspan[0], p[5] is dtspan
    mprob = remake( prob, u0=u0, tspan=tspan, p=p )
    msol  = solve(mprob, LambaEM(), reltol=1e-3,  saveat=dt, callback=cb)
    plot!(msol, lw=2, title="Adaptive SDE", ylim=(0, kmu*1.5) )
  
  
end
 
 
# ---------------

   
@model function fishery_model_turing_sde( M0, kmu, ksd, removed, prob, N=length(M0), ::Type{T} = Float64) where {T}
    # single global model is still to slow to use for parameter estimation
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

    mprob = remake( prob, u0=T[m0[1]*K], tspan=tspan, p=[ r, K ]  )
    msol = solve( mprob, LambaEM(), reltol=1e-3, callback=cb, tstops=survey_time  ) 
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    # process model
    for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        ym = msol.u[j[1]][1] / K
        m0[i] ~ TruncatedNormal( max( ym, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
      end
    end
    
    # observation model
    # @. M0 ~ TruncatedNormal( (m0 *q) + qc, bosd, -5.0, 5.0 )  # M0 in SD units 
    @. M0 ~ TruncatedNormal( (m0 .+ qc) .* q, bosd, 1e-9, 1.2 ) 
end

 
   
@model function fishery_model_turing_incremental_sde( M0, kmu, ksd, removed, prob, N=length(M0), ::Type{T} = Float64) where {T}
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

    t0 = floor(survey_time[1])
    
    # process model
    for i in 2:N
      u0 = T[m0[i-1]*K]
      tsp = (t0+i-1.1, t0+i+0.1)
      mprob = remake( prob, u0=u0, tspan=tsp, p=[ r, K ]  )
      msol = solve( mprob, LambaEM(), reltol=1e-3, callback=cb, tstops=survey_time  ) 
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
    # @. M0 ~ Cauchy( (m0 *q) + qc, bosd  ) 
    @. M0 ~ TruncatedNormal( (m0 .+ qc) .* q, bosd, 1e-9, 1.2 ) 
end

    

#  ----
#  run
   
  # too slow to use right now
  # fmod = fishery_model_turing_sde( M0, kmu, ksd, tspan, prob  )
    
  fmod = fishery_model_turing_incremental_sde( M0, kmu, ksd, tspan, prob  )


  # test run and spin up compilation
  res  =  sample( fmod,  Turing.MH(), 2 )
  res  =  sample( fmod,  Turing.NUTS(0.65), 2   ) 
  
  res  =  sample( fmod,  Turing.PG(30), 10 )

  res  =  sample( fmod,  Turing.MH(), 1000,  thinning=10 , discard_initial=500  )

  # 45 min
  res  =  sample( fmod,  Turing.MH(), MCMCThreads(), 5000, 8, thinning=100, discard_initial=5000 )

  res  =  sample( fmod,  Turing.NUTS(0.65), 1000 , max_depth=12 ) 
 
  
  
  
  # prob = SDEProblem( rs, dprob, LambaEM(),  callback=cb, saveat=dt  ) # tstops=tstops, 
  # res[:,[Symbol("m0[1]") ]]

t0 = floor(survey_time[1])

for u in 1:1000 
  for  i in 1:N
    u0 = [res[u,:K,1] * res[u,Symbol("m0[$i]"),1]]
    p = [res[u,:r,1], res[u,:K,1]]
    tsp = ( t0+i-1.1, t0+i+0.1 )
    msol = solve( 
      remake( prob, u0=u0, tspan=tsp, p=p ), 
      LambaEM(), 
      callback=cb,  
      saveat=dt   ) #
    plot!(msol; alpha=0.05, color=[3 4], xlim=tsp)
  end
end
 
  plot!(; legend=false)

  k0 = [Integer(floor(mean( res[[:"m0[1]"]].value ) *  mean( res[[:K]].value ) ))]
  pm = [mean( res[[:r]].value ), mean( res[[:K]].value ) ]

  msol = solve( remake( prob, u0=k0, tspan=tspan, p=pm ), LambaEM(), saveat=0.1, callback=cb )
  plot!(msol, label="sde-fishing", xlim=tspan)


  msol = solve( remake( prob, u0=k0, tspan=tspan, p=pm ), LambaEM(), saveat=0.1   ) #
  plot!(msol, label="sde-nofishing")

--- check these:
  # back transform M0 to normal scale 
  yhat = ( M0 ./ mean(res[[:"q"]].value) .- mean(res[[:"qc"]].value)) .* mean(res[[:"K"]].value) 
  scatter!(1:N, yhat   ; color=[1 2])
  plot!(1:N, yhat  ; color=[1 2])


  # back transform M0 to normal scale 
  yhat = ( M0  .- mean(res[[:"qc"]].value)) ./ mean(res[[:"q"]].value) .* mean(res[[:"K"]].value) 
  scatter!(survey_time, yhat   ; color=[1 2])
  plot!(survey_time, yhat  ; color=[1 2])
---
  
  w = zeros(N)
  for u in 1:1000  
    for i in 1:N
      w[i] = res[u,:K,1] * res[u,Symbol("m0[$i]"),1]
    end
    plot!(survey_time, w  ;  alpha=0.1, color=[2 2])
  end

  u = zeros(N)
  v = zeros(N)

  for  i in 1:N
    u[i] = mean( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
    v[i] = std( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
  end
  plot!(survey_time, u  ; color=[1 1], lwd=4)
  
  
  # look at predictions:
  si_pred = Vector{Union{Missing, Float64}}(undef, length(M0))
  prob2 = SDEProblem(rs, prob, saveat=dt , callback=cb )
  fmod_pred = fmod( si_pred, kmu, ksd, removed, prob  ) 
 
  predictions = predict(fmod_pred, res)
  y_pred = vec(mean(Array(group(predictions, :M0)); dims = 1));
  
  plot( M0, y_pred )
  sum(abs2, M0 - y_pred) ≤ 0.1



# Generate a MLE estimate.
mle_estimate = optimize(fmod, MLE())

# Generate a MAP estimate.
map_estimate = optimize(fmod, MAP())









plot_autocorr(res; var_names=["r", "K"]);
 

idata = from_mcmcchains( res; library="Turing" )

Plots.plot( survey_time , summarystats(idata.posterior; var_names=["m0"]).mean )
Plots.plot!( survey_time , M0 )


Plots.plot!( survey_time , summarystats(idata.posterior; var_names=["m0"]).mean .* mean( res[[:K]].value ), legend=:bottomright, ylim=(0,55000) )




  # Plot observations.
  msol = is .* K . ./ q .+ qc 

  plot!( fish_time, bm; color=[1 2], linewidth=1)
  scatter!(bm.t, M0'; color=[1 2])


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
 



bm = tzeros(Real, N+M, U)
  
bm[1] ~ TruncatedNormal( b0, bpsd, eps, 1.0)  ;
for i in 2:N 
  o = r * ( 1.0 - bm[i-1] ) ; 
  bm[i] ~ TruncatedNormal(   bm[i-1] * ( 1.0 + o ) - L[i-1]/K, bpsd, eps, 1.0)  ;
end
for i in (N+1):(M+N) 
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
    for i in (ty+1):N 
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
    for i in (ty+1):N 
      ys = ( Y[i,j] / q ) +  qc
      ys ~ TruncatedNormal(  bm[i,j] - L[i,j]/K, bosd, eps, 1.0)  ; #   fall surveys
    end
  end

end  

# deterministic computations: 
F = zeros(sN+M)
B = zeros(N+M)
C = zeros(N+M)

C[1:N] = L[1:N] ./ K
C[(N+1):(M+N)] = er .* bm[(N):(M+N-1)]
C = 1.0 -. C / bm

F =  -log( max.(C, eps) )  ;

   
# -------------------
# parameter estimates for output
MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
BMSY   = exp(K)/2 ; # biomass at MSY
FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY


# recaled estimates

B[1:N] = bm[1:N] *. K - L[1:N] ;
B[(N+1):(M+N)] = (bm[(N+1):(M+N)] - C[(N):(M+N-1)]) *. K ;


 
