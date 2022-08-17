  # testing DDE version -- initial values for size_structured! 
  using DifferentialEquations
  using Interpolations
  using Plots
  
  ks = 100
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* ks
  b=[ 1.0, 0.8 ]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
  d=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
  v=[0.9, 0.9, 0.9, 0.9];  
  tau=1.0 

  survey_time = 1999:2021
  
  external_forcing = ones(length(survey_time),6)  # turns it off
  # external_forcing = rand(length(survey_time),6) # random
  
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


  external_forcing = ones(length(survey_time),6)  # turns it off
  efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
  p = ( b, K, d, v, tau, hsa )   

  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt )
  v = 1; plot( msol.t, reduce(hcat, msol.u)'[:,v], color=[1 1] , alpha=0.75, lw=5 ) 

  plot!( msol, label="dde, no hsa, no fishing", legend=:left )
   

  external_forcing = rand(length(survey_time),6) # random
  efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
  p = ( b, K, d, v, tau, hsa )   

  prob2 = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  # prob2 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )

  msol2 =  solve( prob2,  solver, saveat=dt )
  v = 1; plot!( msol2.t, reduce(hcat, msol2.u)'[:,v], color=[3 3] , alpha=0.75, lw=5 ) 

  plot!( msol2, label="dde, with hsa, no fishing" )

  plot!(; legend=true, xlim=(1999,2021) )
