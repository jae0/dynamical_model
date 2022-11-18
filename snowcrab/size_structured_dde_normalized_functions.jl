using Turing

function size_structured_dde!( du, u, h, p, t)
  # here u, du are actual numbers .. not normalized by K due to use of callbacks

  (b, K, d, d2, v, tau, hsa)  = p

  # u1 = h(p, t-1.0)[2:5]    # no in previous years
  # f8 = h(p, t-8.0)[6]      # no mature fem 8  yrs ago
  # hf = hsa(t, 1:6)

  # this break down seems to speed it up a bit ... not sure why
  br =  h(p, t-8.0)[6]  .* b
  tr =  v .*  h(p, t-1.0)[2:5] 
  
  # d: background mortality + d2: excess mortality due to habitat fluctuations 
  # dr = d .* u  + d2 .* u .* u ./ hf   
  dr = u .* ( d  + d2 .* u ./ hsa(t, 1:6) )  
   
  a = K[2:6] ./ K[1:5]

  du[1] = tr[1] * a[1]            - dr[1]       # note:
  du[2] = tr[2] * a[2]   - tr[1]  - dr[2]
  du[3] = tr[3] * a[3]   - tr[2]  - dr[3]
  du[4] = tr[4] * a[4]   - tr[3]  - dr[4]
  du[5] = br[1] * a[5]   - tr[4]  - dr[5]
  du[6] = br[2]                   - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers

end


function dde_parameters()
    # these are dummy initial values .. just to get things started
    b=[1.7, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    v=[0.65, 0.68, 0.61, 0.79];
    d2 = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    params = ( b, K, d, d2, v, tau, hsa)
    return params
end


 

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM,
    solver=MethodOfSteps(Tsit5()), dt = 0.01, ::Type{T} = Float64) where T

    # basic template

    # biomass process model:
    K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all
 
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )

    qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
    # qc ~ filldist( TruncatedNormal( 0.0, 0.1, -10.0, 10.0), nS )  # uninformative

    model_sd ~ filldist( TruncatedNormal( 0.01, 0.05, 0.01, 1.0 ), nS ) 
 
    # "birth" rate from F(y - 8 to 10)  and for males m5 and femaless
    b ~ filldist( TruncatedNormal(100.0, 1.0, 0.01, 1000.0), 2 )  

    # background mortality
    d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 0.99), nS )

    # excess mortality due to habitat (second order)
    d2 ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

    # transition rates
    v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), 4 )

    # initial conditions
    u0 ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), nS )

    pm = ( b, K, d, d2, v, tau, hsa )
    # @show pm
 
    # process model
    msol = solve(
        remake( prob; u0=u0, h=h, tspan=tspan, p=pm ),
        solver,
        callback=cb,
        # isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) , # this constraint is too much as there are timeseries-lag effects
        isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),  # permit exceeding K
        saveat=dt
    )

    # @show msol.retcode
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood of the data
    for i in 1:nSI
        ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
        for k in 1:nS
            S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
        end
    end

end


@model function size_structured_dde_turing_north( S, kmu, tspan, prob, nT, nS, nM,
  solver=MethodOfSteps(Tsit5()), dt = 0.01, ::Type{T} = Float64) where T

  # biomass process model:
  K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all

  q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )

  qc ~ arraydist([TruncatedNormal( -SminFraction[i], 0.1, -10.0, 10.0) for i in 1:nS])  # informative prior on relative height 
  # qc ~ filldist( TruncatedNormal( 0.0, 0.1, -10.0, 10.0), nS )  # uninformative

  model_sd ~ filldist( TruncatedNormal( 0.0, 0.1, 0.01, 0.3 ), nS ) 

  # "birth" rate from F(y - 8 to 10)  and for males m5 and femaless
  b ~ filldist( TruncatedNormal(10.0, 0.5, 0.01, 100.0), 2 )  

  # background mortality
  d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 0.99), nS )

  # excess mortality due to habitat (second order)
  d2 ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

  # transition rates
  v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), 4 )

  # initial conditions
  u0 ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), nS )

  pm = ( b, K, d, d2, v, tau, hsa )
  # @show pm

  # process model
  msol = solve(
      remake( prob; u0=u0, h=h, tspan=tspan, p=pm ),
      solver,
      callback=cb,
      # isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) , # this constraint is too much as there are timeseries-lag effects
      isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),  # permit exceeding K
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
      end
  end

end


@model function size_structured_dde_turing_south( S, kmu, tspan, prob, nT, nS, nM,
  solver=MethodOfSteps(Tsit5()), dt=0.01, ::Type{T} = Float64) where {T}

  # biomass process model:
  K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all

  q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )

  qc ~ arraydist([TruncatedNormal( -SminFraction[i], 0.1, -10.0, 10.0) for i in 1:nS])  # informative prior on relative height 
  # qc ~ filldist( TruncatedNormal( 0.0, 0.1, -10.0, 10.0), nS )  # uninformative

  model_sd ~ filldist( TruncatedNormal( 0.0, 0.1, 0.01, 0.3 ), nS ) 

  # "birth" rate from F(y - 8 to 10)  and for males m5 and femaless
  b ~ filldist( TruncatedNormal(10.0, 0.5, 0.01, 100.0), 2 )  

  # background mortality
  d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 0.99), nS )

  # excess mortality due to habitat (second order)
  d2 ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

  # transition rates
  v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), 4 )

  # initial conditions
  u0 ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), nS )

  pm = ( b, K, d, d2, v, tau, hsa )
  # @show pm

  # process model
  msol = solve(
      remake( prob; u0=u0, h=h, tspan=tspan, p=pm ),
      solver,
      callback=cb,
      # isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) , # this constraint is too much as there are timeseries-lag effects
      isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),  # permit exceeding K
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
      end
  end

end


@model function size_structured_dde_turing_4x( S, kmu, tspan, prob, nT, nS, nM,
  solver=MethodOfSteps(Tsit5()), dt = 0.01, ::Type{T} = Float64) where T

  # biomass process model:
  K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all

  q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )

  qc ~ arraydist([TruncatedNormal( -SminFraction[i], 0.1, -10.0, 10.0) for i in 1:nS])  # informative prior on relative height 
  # qc ~ filldist( TruncatedNormal( 0.0, 0.1, -10.0, 10.0), nS )  # uninformative

  model_sd ~ filldist( TruncatedNormal( 0.0, 0.1, 0.01, 0.3 ), nS ) 

  # "birth" rate from F(y - 8 to 10)  and for males m5 and femaless
  b ~ filldist( TruncatedNormal(10.0, 0.5, 0.01, 100.0), 2 )  

  # background mortality
  d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 0.99), nS )

  # excess mortality due to habitat (second order)
  d2 ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

  # transition rates
  v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), 4 )

  # initial conditions
  u0 ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), nS )

  pm = ( b, K, d, d2, v, tau, hsa )
  # @show pm

  # process model
  msol = solve(
      remake( prob; u0=u0, h=h, tspan=tspan, p=pm ),
      solver,
      callback=cb,
      # isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) , # this constraint is too much as there are timeseries-lag effects
      isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),  # permit exceeding K
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
      end
  end

end


# ------------------------------


function fishery_model_test( test=("basic", "random_external_forcing", "fishing", "nofishing") )

  ## test DifferentialEquations DDE model 
  gr()
  theme(:default)
 
  if any( occursin.( r"basic", test )  )

    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    b=[0.8, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* kmu;
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    v=[0.65, 0.68, 0.61, 0.79];
    d2 =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
    
    u0 = [ 0.65, 0.6, 0.52, 0.62, 0.58, 0.32 ]  ; 
    tau=[1.0] 
    params = ( b, K, d, d2, v, tau, hsa)
    prob1 = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]
    out = msol1 =  solve( prob1,  solver, callback=cb, saveat=dt )
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
  
    p = ( b, K, d, d2, v, tau, hsa )   
    tspan = (1990.0, 2050.0)
    nS = length(u0)  # n components
    
    # history function for time before t0:: 0.5 default
    # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS)  .*0.5
    tau = 1.0  # delay
    
    solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
    prob2 = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=tau )
    out = msol2 =  solve( prob2,  solver, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol2 ; title="basic" )
  
  elseif  any( occursin.( r"nofishing", test ) )

    ks = 1.0e9
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    b=[0.6, 0.4]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    d2=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6] 
    tau=[1.0] 

    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
   
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    external_forcing = ones(length(survey_time),6)  # turn off external forcing  (wrt viable habitat)
    efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
    p = ( b, K, d, d2, v, tau, hsa )   
  
    prob3 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    out = msol3 =  solve( prob3,  solver, saveat=dt  ) #, isoutofdomain=(y,p,t)->any(x->(x<0)|(x>1), y) )
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
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
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
   
    p = ( b, K, d, d2, v, tau, hsa)   
  
    prob4 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    # prob4 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
      
    out = msol4 =  solve( prob4,  solver, saveat=dt )
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

  

function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=100 )

  nchains = size(res)[3]
  nsims = size(res)[1]

  nZ = nchains*nsims
  nI = Int( min( nZ , n_sample ) )
 
  md = zeros(nM, nS, nZ, 2)  # number normalized
  mn = zeros(nM, nS, nZ, 2)  # numbers
  mb = mn[:,1,:,:]  # biomass of first class
  z = 0

  for j in 1:nsims  # nsims
  for l in 1:nchains #nchains
    z += 1
    
    b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    d2=[ res[j, Symbol("d2[$k]"), l] for k in 1:nS]

    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

    pm = ( b, K, d, d2, v, tau, hsa )

    prb = remake( prob; u0=u0 , h=h, tspan=tspan, p=pm )

    msol = solve( prb, solver, callback=cb, saveat=dt )
    msol2 = solve( prb, solver, saveat=dt ) # no call backs

    
    for i in 1:nM
      
        ii = findall(x->x==prediction_time[i], msol.t) 
        jj = findall(x->x==prediction_time[i], msol2.t)
        
        if length(ii) > 0 && length(jj) > 0 
          ii = ii[1]
          jj = jj[1]
          sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t[ii])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
          sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t[jj]) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
          md[i,:,z,1] = msol.u[ii]   # with fishing
          md[i,:,z,2] = msol2.u[jj]  # no fishing
          mn[i,:,z,1] = msol.u[ii]   .* K  # with fishing
          mn[i,:,z,2] = msol2.u[jj]   .* K # no fishing
          mb[i,z,1] = mn[i,1,z,1]  .* sf
          mb[i,z,2] = mn[i,1,z,2]  .* sf2
        end

    end

  end
  end

  # extract sims (with fishing)
  g = mb[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
  # size(g)

  # plot biomass
  gr()
  pl = plot()
  pl = plot!(pl, prediction_time, g[:,sample(1:nZ, nI)];  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )

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
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )

  return (md, mn, mb, pl)

end


# -----------

function abundance_from_index( aindex, res, k )
  u = size(res)
  v = aindex .* ones(length(aindex), u[1]) # expand to matrix
  yhat = ( v' .- res[:,Symbol("qc[$k]"),:]  )  ./  res[:,Symbol("q[$k]"),:] .* res[:,Symbol("K[$k]"),:]  
  return yhat'
end

# -----------


function fishery_model_predictions_trace( res; n_sample=10, plot_k=1, alpha=0.02, plot_only_fishing=true )

    nchains = size(res)[3]
    nsims = size(res)[1]

    nZ = nchains*nsims
    nI = Int( min( nZ , n_sample ) )
    
    jj = sample( 1:nsims, nI )  # to plot
    ll = sample( 1:nchains, nI )  # to plot
   
    out = Vector{Vector{Float64}}()
    out2 = Vector{Vector{Float64}}()

    gr()
    theme(:default)
    pl =plot()

    for z in 1:nI  # nsims
        j = jj[z]
        l = ll[z]

        b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
        K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
        d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
        v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
        d2=[ res[j, Symbol("d2[$k]"), l] for k in 1:nS]

        u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

        pm = ( b, K, d, d2, v, tau, hsa )

        prb = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

        if plot_k==1
            # do fishing and nonfishing

            msol = solve( prb, solver, callback=cb, saveat=dt )
            msol2 = solve( prb, solver, saveat=dt ) # no call backs

            sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t) ./ 1000.0 ./ 1000.0 :  scale_factor
            sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt

            yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k]) .* K[plot_k] .* sf2
            yval = vec( reduce(hcat, msol.u)'[:,plot_k] ) .* K[plot_k] .* sf

            pl = plot!( pl, msol.t, yval, alpha=alpha, lw=1, color=:orange )
            if !plot_only_fishing
              pl = plot!( pl, msol2.t, yval2, alpha=alpha*2.0, lw=1, color=:lime )
            end

            push!(out, yval)
            push!(out2, yval2)

        else
            msol2 = solve( prb, solver, saveat=dt ) # no call backs
            yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k]) .* K[plot_k]
            pl = plot!( pl, msol2.t, yval2, alpha=alpha, lw=1, color=:lime )

            push!(out2, yval2)
        end
 
    end

    pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
    # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
    pl =  plot!(pl; legend=false )

    return (out, out2, pl)
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



function fishery_model_harvest_control_rule(res, yrs; FM=FM, fb=fb, n_sample=500 )

  fmsy = nothing

  pl = plot()

  # mean weight by year
  sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor

  # sample and plot posterior K
  K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 

  o = sample(K, n_sample)
  pl = vline!(pl, o;  alpha=0.05, color=:limegreen )
  pl = vline!(pl, o./2;  alpha=0.05, color=:darkkhaki )
  pl = vline!(pl, o./4;  alpha=0.05, color=:darkred )

  pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
  pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

  pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
  pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )

  pl = vline!(pl, [mean(o)/4.0];  alpha=0.6, color=:darkred, lw=5 )
  pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )

  nt = length(survey_time)
  colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]

  # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)

  fb_mean = mean(fb, dims=2)
  fm_mean = mean(FM, dims=2)

  # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
  pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)

  pl = scatter!(pl,  fb_mean, fm_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
    series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
  pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
  
  ub = max( quantile(o, 0.95), maximum( fb_mean ) ) * 1.05
  pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
  # TODO # add predictions ???

  return(K, fb_mean, fm_mean, fmsy, pl)
end




function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end


function removals_aggregate( removed, fish_year )
  landings_aggregated = DataFrame( yr=floor.(fish_year), rem = removed );
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


function plots_diagnostic( res, vn="K" ) 
  gr()
  pl = plot()
  pl = density!(pl, res[ Symbol(vn) ])
  return pl
  # vn = "b[1]"; density!(res[ Symbol(vn) ])
  # vn = "b[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[1]"; density!(res[ Symbol(vn) ])
  # vn = "d[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[3]"; density!(res[ Symbol(vn) ])
  # vn = "d[4]"; density!(res[ Symbol(vn) ])
  # vn = "d[5]"; density!(res[ Symbol(vn) ])
  # vn = "d[6]"; density!(res[ Symbol(vn) ])
  # vn = "K[1]"; density!(res[ Symbol(vn) ])
  # vn = "K[2]"; density!(res[ Symbol(vn) ])
  # vn = "K[3]"; density!(res[ Symbol(vn) ])
  # vn = "K[4]"; density!(res[ Symbol(vn) ])
  # vn = "K[5]"; density!(res[ Symbol(vn) ])
  # vn = "K[6]"; density!(res[ Symbol(vn) ])
end



# ----------


function fishery_model_inference( fmod; rejection_rate=0.65, n_adapts=1000, n_samples=1000, n_chains=1, max_depth=7, init_ϵ=0.05, 
  turing_sampler = Turing.NUTS(n_adapts, rejection_rate; max_depth=max_depth, init_ϵ=init_ϵ), seed=1  )
   
  Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
  
  Random.seed!(seed)
  
  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs;; 29hrs 100/100 cfasouth
  #   # n_chains = Threads.nthreads()
  # turing_sampler = Turing.NUTS(n_adapts, rejection_rate; max_depth=max_depth, init_ϵ=init_ϵ)  ;# stepsize based upon previous experience
  
  res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains )
  # if on windows and threads are not working, use single processor mode:
  # res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

  showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

  return res
end



# -----------


function fishery_model_mortality( removed, fb; n_sample=100 )    

  removed_annual_kt = removals_aggregate( removed, fish_year )
  Fkt = removed_annual_kt[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt

  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  # FM[ FM .< eps(0.0)] .= zero(eltype(FM))

  o = mean( FM, dims=2)
  o[isnan.(o)] .= zero(eltype(FM))
 
  ub = maximum(o) * 1.1
  nchains = size(res)[3]
  nsims = size(res)[1]
  nZ = nchains*nsims
  nI = Int( min( nZ , n_sample ) )

  pl = plot()
  pl = plot!(pl, survey_time, FM[:,sample(1:nZ, nI)] ;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, survey_time, o ;  alpha=0.8, color=:slateblue, lw=4)
  pl = plot!(pl, xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  pl = plot!(pl, ylim=(0, ub ) )
  pl = plot!(pl ; legend=false )
  return ( Fkt, FR, FM, pl )
end

