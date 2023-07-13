using Turing
  

function size_structured_dde!( du, u, h, p, t )
  # here u, du are actual numbers .. not normalized by K due to use of callbacks

  # b, K, d, d2, v, tau, hsa  = p
  b, K, d, d2, v  = p  # hsa is global 
 
  @inbounds begin

      # this break down seems to speed it up a bit ... not sure why
      br =  b .* h(p, t-8.0)[6]  
      tr =  v .* h(p, t-1.0)[2:5]
      dh = d2 .* u ./ hsa(t, 1:6)
      dr = d .* u  .+  dh .* u 
      
      du[1] = tr[1] * K[2] / K[1]            - dr[1]       # note:
      du[2] = tr[2] * K[3] / K[2]   - tr[1]  - dr[2]
      du[3] = tr[3] * K[4] / K[3]   - tr[2]  - dr[3]
      du[4] = tr[4] * K[5] / K[4]   - tr[3]  - dr[4]
      du[5] = br[1] * K[6] / K[5]   - tr[4]  - dr[5]
      du[6] = br[2]                          - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers
    
    end

end


function dde_parameters()
    # these are dummy initial values .. just to get things started
    b=[1.7, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    d2 = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    v=[0.65, 0.68, 0.61, 0.79];
    # params = ( b, K, d, d2, v, tau, hsa)
    params = ( b, K, d, d2, v)
    return params
end


function firstindexin(a::AbstractArray, b::AbstractArray)
  bdict = Dict{eltype(b), Int}()
  for i=length(b):-1:1
      bdict[b[i]] = i
  end
  [get(bdict, i, 0) for i in a]
end
 

function β( mode, conc )
  # alternate parameterization of beta distribution 
  # conc = α + β     https://en.wikipedia.org/wiki/Beta_distribution
  beta1 = mode *( conc - 2  ) + 1.0
  beta2 = (1.0 - mode) * ( conc - 2  ) + 1.0
  Beta( beta1, beta2 ) 
end 

 
@model function size_structured_dde_turing( ; PM, solver_params )
 
  K ~ filldist( LogNormal( PM.logkmu[1], PM.logkmu[2]), PM.nS )  # kmu already on log scale # is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( PM.q[1], PM.q[2] ), PM.nS )
  qc ~ arraydist([Normal( PM.qc[1][i], PM.qc[2]) for i in 1:PM.nS])  # informative prior on relative height 
  model_sd ~  arraydist(LogNormal.( PM.logScv[1], PM.logScv[2]) ) # working: β(0.1, 10.0);  plot(x->pdf(β(0.3,12), x), xlim=(0,1)) # uniform 
  b ~   filldist( LogNormal( PM.b[1],  PM.b[2] ), PM.nB )   # centered on 1; plot(x->pdf(LogNormal(log(10), 1.0), x), xlim=(0,10)) # mode of 5
  d ~   filldist( LogNormal( PM.d[1],  PM.d[2] ), PM.nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( PM.d2[1], PM.d2[2]), PM.nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 
  v ~   filldist( LogNormal( PM.v[1],  PM.v[2] ), PM.nG ) # transsition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  
  u0 ~  filldist( Beta(2, 2), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # process model
  prob_new = remake( solver_params.prob; u0=u0, h=solver_params.h, 
    tspan=solver_params.tspan, p=( b, K, d, d2, v ) )

  msol = solve(
      prob_new,
      solver_params.solver, 
      callback=solver_params.cb,
      abstol=solver_params.abstol, 
      reltol=solver_params.reltol, 
      # tstops=solver_params.saveat,  
      saveat=solver_params.dt
      # ,
      # saveat=solver_params.saveat
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  ii = indexin(PM.Stime,  msol.t)
  ii = ii[ .!(isnothing.(ii)) ]

  if length(ii) != PM.nSI
    Turing.@addlogprob! -Inf
    return nothing
  end

  
  # likelihood of the data
  A = view( Array(msol), :, ii) .* q .+ qc
  sigma = view(  model_sd, PM.Si, : )
  
  @. PM.data ~ Normal( A', sigma )  
  # @. PM.data ~ Normal( A', model_sd[PM.Si, :] )   
  # PM.datavector ~ MvNormal( vec(A'), Diagonal(vec(sigma).^2.0)  )  # no real speed gain using MVN .. JC: Feb 2023

end
 
# ------------------------------
  

function fishery_model_test( test=("basic", "random_external_forcing", "fishing", "nofishing") )

  ## test DifferentialEquations DDE model 
  gr()
  theme(:default)
  solver_test =  MethodOfSteps(Tsit5())

  if any( occursin.( r"basic", test )  )

    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    
    b=[0.8, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* kmu;
    v=[0.65, 0.68, 0.61, 0.79];
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    d2 =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
    
    u0 = [ 0.65, 0.6, 0.52, 0.62, 0.58, 0.32 ]  ; 
    tau=[1.0] 
    # params = ( b, K, d, d2, v, tau, hsa )
    params = ( b, K, d, d2, v  )
 
    prob1 = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]

    out = msol1 =  solve( prob1, solver_test, callback=cb, saveat=dt )
    pl = plot()
    pl = plot( pl, msol1, ; legend=false, xlim=(1999,2021), title="basic" )

  elseif  any( occursin.( r"random_external_forcing", test )  )

    # testing DDE version -- initial values for size_structured_dde! 
    # basic run
    ks = 1000
    u0 = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ] ;
    b=[ 0.6, 0.6 ]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    d=[0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
    v=[0.9, 0.9, 0.9, 0.9];  
    d2=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
    tau=[1.0] 
  
    survey_time = 1999:2021
    
    external_forcing = ones(length(survey_time),6)  # turns it off
    # external_forcing = rand(length(survey_time),6) # random
    
    efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc, 1999:2021, 1:6 )
    
    # p = ( b, K, d, d2, v, tau, hsa )   
    p = ( b, K, d, d2, v )   
    
    tspan = (1990.0, 2050.0)
    nS = length(u0)  # n components
    
    # history function for time before t0:: 0.5 default
    # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS)  .*0.5
    tau = [1.0]  # delay
    
    prob2 = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=tau )
    out = msol2 =  solve( prob2,  solver_test, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol2 ; title="basic" )
  
  elseif  any( occursin.( r"nofishing", test ) )

    ks = 1.0e9
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    b=[0.6, 0.4]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6] 
    tau=[1.0] 

    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
   
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    external_forcing = ones(length(survey_time),6)  # turn off external forcing  (wrt viable habitat)
    efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
    
    # p = ( b, K, d, d2, v, tau, hsa )   
    p = ( b, K, d, d2, v )   
    
    prob3 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    
    out = msol3 =  solve( prob3,  solver_test, saveat=dt  ) #, isoutofdomain=(y,p,t)->any(x->(x<0)|(x>1), y) )
    pl = plot()

    pl = plot!( pl, msol3, title="dde, with random hsa, with fishing") 

    i = 1; # fb
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_1" ) 

    i = 6; # fem
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_6" ) 
  
   
  elseif  any( occursin.( r"fishing", test )  | occursin( r"nofishing", test ) )

    ks = 1.0e9
    K = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    b=[0.8, 0.5]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2=[0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    tau=[1.0] 

    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks
   
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    # add random viable habitat forcing
    external_forcing = rand(length(survey_time),6) # random
    efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
   
    #p = ( b, K, d, d2, v, tau, hsa)   
    p = ( b, K, d, d2, v )   
   
    prob4 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    # prob4 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
      
    out = msol4 =  solve( prob4,  solver_test, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol4, title="dde, with random hsa, no fishing" )

    i = 1; 
    pl = plot!( pl, msol4.t, reduce(hcat, msol4.u)'[:,i], color=[3 3] , alpha=0.75, lw=5 ) 
    
    i = 6; 
    pl = plot!( pl, msol4.t, reduce(hcat, msol4.u)'[:,i], color=[1 1] , alpha=0.75, lw=5 ) 
    
  end

  pl = plot!( pl; legend=true, xlim=(1990,2021) )

  return (msol, pl) 

end
 
 

# -----------

function abundance_from_index( aindex, res, k )
  u = size(res)
  v = aindex .* ones(length(aindex), u[1]) # expand to matrix
  yhat = ( v' .- res[:,Symbol("qc[$k]"),:]  )  ./  res[:,Symbol("q[$k]"),:] .* res[:,Symbol("K[$k]"),:]  
  return yhat'
end


# -----------

  


function fishery_model_predictions_timeseries( num; prediction_time, plot_k )
  gk = num[:,plot_k,:,1]
  pl = plot()
  pl = plot!(pl, prediction_time, gk;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  yhat = ( S[:,plot_k] .- mean(res[:,Symbol("qc[$plot_k]"),:]  ) ) ./ mean(res[:,Symbol("q[$plot_k]"),:]) .* mean(res[:,Symbol("K[$plot_k]"),:]  )

  pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )
  return (gk, pl)
end

# -----------





function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end


function removals_aggregate( removed, fish_year )
  landings_aggregated = DataFrame( yr=floor.(fish_year), remNum = removed );
  landings_aggregated.rem = landings_aggregated.remNum .* mw(fish_time) 
  out = combine(groupby(landings_aggregated,:yr),[:rem ] .=> sum )
  stimes = DataFrame( yr=floor.(survey_time) )
  out = leftjoin(stimes, out, on=:yr)
  sort!(out, :yr)
  oo = findall(x->ismissing(x), out[:,:rem_sum])
  if length(oo) > 0
    out[ oo[1], :rem_sum ] = 0.0
  end
  return(out)
end


function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end


function plots_diagnostic( res; vn="K", i=-1, toplot="" ) 
  gr()
  pl = plot()
  if toplot == "" 
    if i != -1 
      toplot = Symbol(vn,"[$i]")
    else 
      toplot = Symbol(vn)
    end
  end
  pl = density!(pl, res[ toplot ])
  return pl
 end


# ----------

  

function expand_grid(; kws...)
  names, vals = keys(kws), values(kws)
  return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end
 


# -------------------


function fishery_model_predictions( res; prediction_time=prediction_time, solver_params=solver_params, PM=PM, 
  n_sample=-1, lower_bound=0.0, override_negative_solution=false, ntries_mult=10 )

  nchains = size(res)[3]
  nsims = size(res)[1]
  
  if n_sample == -1
    # do all
    n_sample = nchains * nsims
  end

  md = zeros(nM, nS, n_sample, 2)  # number normalized
  mn = zeros(nM, nS, n_sample, 2)  # numbers
  mb = mn[:,1,:,:]  # biomass of first class
 
  trace_time = collect( solver_params.tspan[1]:solver_params.dt:solver_params.tspan[2] )
  ntt = length(trace_time)
   
  trace_bio = zeros(3, ntt, n_sample ) 
  trace = zeros( nS, ntt, n_sample )

  ntries = 0
  z = 0

  while z <= n_sample 
    ntries += 1
    ntries > ntries_mult * n_sample && break 
    z >= n_sample && break

    j = rand(1:nsims)  # nsims
    l = rand(1:nchains) # nchains

    b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
    d2=[ res[j, Symbol("d2[$k]"), l] for k in 1:nS]
    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

    prb = remake( solver_params.prob; u0=u0 , h=solver_params.h, tspan=solver_params.tspan, p=( b, K, d, d2, v ) )
    msol1 = solve( prb, solver_params.solver, callback=solver_params.cb, 
      saveat=solver_params.dt ,
      isoutofdomain=(y,p,t)->any(x -> x<lower_bound, y)  # permit exceeding K
    )

    msol0 = solve( prb, solver_params.solver,  
      saveat=solver_params.dt,
      isoutofdomain=(y,p,t)->any(x -> x<lower_bound, y)  # permit exceeding K
    )

    if msol1.retcode == :Success && msol0.retcode == :Success
     
      ii0 = firstindexin( prediction_time, msol0.t )
      ii1 = firstindexin( prediction_time, msol1.t )
      
      jj0 = firstindexin( trace_time, msol0.t)
      jj1 = firstindexin( trace_time, msol1.t)
     
      if length(ii0) != PM.nM | length(ii1) != PM.nM
        continue
      end
      
      if length(jj0) != ntt | length(jj1) != ntt
        continue
      end
      
      z += 1
 
      MS0 = Array(msol0) 
      MS1 = Array(msol1) 

      if override_negative_solution
        MS0[ findall(x -> x<0.0, skipmissing(MS0) ) ] .= 0.0
        MS1[ findall(x -> x<0.0, skipmissing(MS1) ) ] .= 0.0
      end

      sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol1.t[ii1])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
      
      md[:,:,z,1] = MS1[:,ii1]'  # with fishing
      md[:,:,z,2] = MS0[:,ii0]'  # witout fishing
      mn[:,:,z,1] = MS1[:,ii1]'  .* K' # with fishing scaled to K
      mn[:,:,z,2] = MS0[:,ii0]'  .* K' # without fishing scaled to K
      mb[:,z,1] = mn[:,1,z,1]  .* sf  # biomass of state var 1 
      mb[:,z,2] = mn[:,1,z,2]  .* sf

      # traces for plotting, etc 
      sft  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(trace_time)  ./ 1000.0 ./ 1000.0 :  scale_factor
      trace_bio[1,:,z] = MS1[1,jj1] .* K[1] .* sft  # with 
      trace_bio[2,:,z] = MS0[1,jj0] .* K[1] .* sft  # without
      trace_bio[3,:,z] = ( trace_bio[2,:,z] .- trace_bio[1,:,z] ) ./ trace_bio[2,:,z] 
 
      trace[:,:,z] = MS1[:,jj1] .* K
    end # if
  end  # while
 
  if z < n_sample 
    @warn  "Insufficient number of solutions" 
  end
 
  sols= findall( x -> x!=0, vec( sum(mb, dims=(1,3,4) )) )
  if length(sols) > 1 
    md = md[:,:,sols,:]
    mn = mn[:,:,sols,:]
    mb = mb[:,sols,:,:]
  end
 
  return (md, mn, mb, trace, trace_bio, trace_time )

end



# -----------


function fishery_model_mortality( ; removed=removed, bio=bio, survey_time=survey_time, fish_year=fish_year )   
  fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete 
  removed_annual_kt = removals_aggregate( removed, fish_year )
  Fkt = removed_annual_kt[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  # FM[ FM .< eps(0.0)] .= zero(eltype(FM))
  return ( Fkt, FR, FM  )
end

# -----------


function fishery_model_plot(; toplot=("fishing", "nofishing", "survey"), n_sample=200,
  res=res, bio=bio, num=num, trace=trace, trace_bio=trace_bio, FM=FM, 
  S=S, si=1, scale_factor=scale_factor, 
  prediction_time=prediction_time, survey_time=survey_time, yrs=yrs, 
  alphav=0.05, 
  pl= Plots.plot(), 
  time_range=(floor(minimum(survey_time))-0.25, ceil(maximum(survey_time))+0.25 ),
  time_range_predictions=(floor(minimum(survey_time))-1.0, ceil(maximum(prediction_time)) )
)

  trace_time = collect( solver_params.tspan[1]:solver_params.dt:solver_params.tspan[2] )
   
  nsims = size(bio)[2]
  ss = rand(1:nsims, n_sample)  # sample index

  if any(isequal.("nofishing", toplot))  
    g = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:lime)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("fishing", toplot))  
    g = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("footprint", toplot))  
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    pl = plot!(pl, prediction_time, g[:,ss] ;  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
    
  end

  if any(isequal.("survey", toplot))  
    # back transform S to normal scale .. do sims too (TODO)
    k = 1
    qcm = mean(res[:,Symbol("qc[$k]"),:])
    qm  = mean(res[:,Symbol("q[$k]"),:])
    Km = mean(res[:,Symbol("K[$k]"),:]  )
    yhat = ( S[:,k] .- qcm  ) ./ qm .* Km   # abundance_from_index S[:,1]  
    if nameof(typeof(mw)) == :ScaledInterpolation
      yhat = yhat .* mw(yrs) ./ 1000.0  ./ 1000.0
    else
      yhat = yhat .* scale_factor
    end
    pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, yhat, markersize=4, color=:darkgray)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )

  end
  

  if any(isequal.("number", toplot))  

    gk = num[:,si,:,1]
    pl = plot!(pl, prediction_time, gk[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

    # back transform S to normal scale .. do sims too (TODO)
    yhat = ( S[:,si] .- mean(res[:,Symbol("qc[$si]"),:]  ) ) ./ mean(res[:,Symbol("q[$si]"),:]) .* mean(res[:,Symbol("K[$si]"),:]  )

    pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("number_trace", toplot))  

    gk = trace_bio[si,:,:]
    pl = plot!(pl, trace_time, gk[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, trace_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

    # back transform S to normal scale .. do sims too (TODO)
    yhat = ( S[:,si] .- mean(res[:,Symbol("qc[$si]"),:]  ) ) ./ mean(res[:,Symbol("q[$si]"),:]) .* mean(res[:,Symbol("K[$si]"),:]  )

    pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("trace", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl; legend=false )
    pl = plot!( pl; xlim=time_range )
  end

  if any(isequal.("trace_projections", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl, trace_time, mean(trace_bio[2,:,ss], dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!( pl, trace_time, mean(trace_bio[1,:,ss], dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!( pl; legend=false )
    pl = plot!( pl; xlim=time_range_predictions )
    pl = vline!( pl, year_assessment .+ [1,2],  alpha=0.5, color=:lightslategray, line=:dash, lw=1.5 )
  end

  if any(isequal.("trace_nofishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2,:,ss], alpha=alphav, lw=1, color=:lime )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_fishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[1,:,ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end
 
  if any(isequal.("trace_footprint", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3,:,ss], alpha=alphav, lw=1, color=:lightslateblue )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_footprint_projections", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3,:,ss], alpha=alphav, lw=1, color=:lightslateblue )
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=(floor(minimum(survey_time))-1.0, ceil(maximum(prediction_time)) ) )
    pl = vline!( pl, year_assessment .+ [1,2],  alpha=0.5, color=:lightslategray, line=:dash, lw=1.5 )
    pl = plot!(pl; xlim=time_range_predictions )
  end

  if any(isequal.("fishing_mortality", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, survey_time, FM[:,ss] ;  alpha=0.02, color=:lightslateblue)
    pl = plot!(pl, survey_time, FMmean ;  alpha=0.8, color=:slateblue, lw=4)
    pl = plot!(pl, ylim=(0, ub ) )
    pl = plot!(pl ; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("fishing_mortality_vs_footprint", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    g = g[1:length(survey_time),:]
    pl = scatter!(pl, FM[:,ss], g[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = scatter!(pl, FMmean, mean(g, dims=2);  
      alpha=0.8, color=:darkslateblue, lw=4, markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4))
    pl = plot!(pl ; legend=false )
  end


  if any(isequal.("harvest_control_rule_footprint", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
  
    fb_mean = mean(fb, dims=2)

    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
  
    g = g[1:length(survey_time),:]
  
    g_mean = mean(g, dims=2)
  
    # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
    pl = plot!(pl, fb_mean, g_mean ;  alpha=0.8, color=:slateblue, lw=3)
  
    pl = scatter!(pl,  fb_mean, g_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
    pl = scatter!(pl,  [fb_mean[nt]], [g_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
    
    ub = max( quantile(K, 0.95), maximum( fb_mean ) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(g_mean ) * 1.05  ) )
 
  end
   

  if any(isequal.("harvest_control_rule", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
  
    fb_mean = mean(fb, dims=2)
    fm_mean = mean(FM, dims=2)
  
    # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
    pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)
  
    pl = scatter!(pl,  fb_mean, fm_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
    pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
    
    ub = max( quantile(K, 0.95), maximum( fb_mean ) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
 
  end
   
  return(pl)

end
 

function fishing_pattern_from_data(  fish_time, removed, ny=5 ) 
  # choose last n years of fishing and model attack rate 
  # scale to 1 then rescale to % of catch by time

  f = DataFrame( fish_time=fish_time, removals=removed )  # number
  f.yr = floor.( f.fish_time)
  f = f[(f[!,:yr] .> maximum(f[!,:yr]) - ny ),:]
  
  annual = DataFrame(
    fsum = [sum(x[!,:removals]) for x in groupby(f, :yr)],
    yr =   [ x[1,:yr] for x in groupby(f, :yr)]
  )

  f = innerjoin( f, annual, on=:yr)
  f.sy = f.fish_time - f.yr
  dtt = 12
  f.sy = round.( floor.( f.sy * dtt ) / dtt, digits=3) 
  f.rem = f.removals ./ f.fsum
 
  allowmissing!(f)
  f[findall(x->x>0.99, f.rem),:rem] .= missing

  ff = DataFrame(
    sy = [ x[1,:sy] for x in groupby(f, :sy)],
    fa = [mean(x[!,:rem ]) for x in groupby(f, :sy)]
  )
  sort!(ff, [order(:sy)])
  ff.fa = ff.fa / sum(ff.fa)

  gg = DataFrame( mon=0:1:12 )
  gg.sy = round.( floor.( gg.mon / 12 * dtt ) / dtt, digits=3) 
  
  fs = outerjoin( gg, ff, on=:sy)
  fs[findall(x->ismissing(x), fs.fa),:fa] .= 0
  sort!(fs, [order(:sy)])
  
  ifs = extrapolate( interpolate( fs.fa, (BSpline(Linear())) ), Interpolations.Flat() )
  fishing_pattern_seasonal_interpolation_function = Interpolations.scale(ifs, 0:1/12:1 )

  return fishing_pattern_seasonal_interpolation_function
end


function project_with_constant_catch( res; solver_params=solver_params, PM=PM, Catch=0, ny_fishing_pattern=5 )
  
  # forward project assuming constant fishing pattern
  # callbacks for external perturbations to the system (deterministic fishing without error)
   
  sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
  # sample and plot posterior K
  K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  ER =  Catch / mean(K) 

  exploitationrate = exp(ER)-1.0  # relative to K

  fishing_pattern_seasonal = fishing_pattern_from_data(fish_time, removed, ny_fishing_pattern ) # fraction of annual total .. fishing_pattern(0.1) gives fraction captured on average by 0.1 * 365 days 
  
  condition_fp = function(u, t, integrator )
    t in fish_time_project
  end

  function affect_fishing_project!(integrator)
    k = integrator.t - floor(integrator.t)
    integrator.u[1] -=  exploitationrate / hsa(integrator.t,1) * fishing_pattern_seasonal(k)  # p[2] ==K divide by K[1]  .. keep unscaled to estimate magnitude of other components
  end
 
  sp = deepcopy(solver_params)
  sp = @set sp.cb =  CallbackSet(
    PresetTimeCallback( fish_time, affect_fishing! ),
    DiscreteCallback( condition_fp, affect_fishing_project!, save_positions=(true, true) )
  );
  
  # n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
  # override_negative_solution=true makes it more permissive .. but tuncates at 0, so be careful
  m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(
    res, solver_params=sp, PM=PM , 
    lower_bound=-0.05, override_negative_solution=true, ntries_mult=10 
  )

  return m, num, bio, trace, trace_bio, trace_time
end


 