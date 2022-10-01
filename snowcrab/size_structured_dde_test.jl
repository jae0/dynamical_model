  # testing DDE version -- initial values for size_structured_dde! 
  using DifferentialEquations
  using Interpolations
  using Plots

  gr()
  theme(:default)


  ks = 1000
  u0 = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ]  .* ks;
  b=[ 0.75, 0.6 ]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
  d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.3];
  dh=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
  v=[0.9, 0.9, 0.9, 0.9];  
  tau=1.0 

  survey_time = 1999:2021
  
  external_forcing = ones(length(survey_time),6)  # turns it off
  # external_forcing = rand(length(survey_time),6) # random
  
  efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc, 1999:2021, 1:6 )

  p = ( b, K, d, dh, v, tau, hsa )   
  tspan = (1990.0, 2050.0)
  nS = length(u0)  # n components
  # history function for time before t0:: 0.5 default
  # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks
  tau = 1.0  # delay
  lags = [tau]
  solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
  prob = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=0.001 )

  gr()
  plot()
  plot!( msol ; legend=true)
  
  # ---------------
  #test runs
  ks = 1.0e9
  u0 = [ 0.2, 0.5, 0.3, 0.2, 0.1, 0.1 ]  .* ks;;
  b=[ 0.5, 0.4 ]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
  d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
  dh=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

  v=[0.9, 0.9, 0.9, 0.9];  

  h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks

  tau=1.0 
  tspan = (1990.0, 2025.0)

  dt = 0.001
  nS = length(u0)  # n components
 
  external_forcing = ones(length(survey_time),6)  # turns it off
  efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
  p = ( b, K, d, dh, v, tau, hsa )   

  prob = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=lags )
  msol =  solve( prob,  solver, saveat=dt  ) #, isoutofdomain=(y,p,t)->any(x->(x<0)|(x>1), y) )

  plot(0)
  i = 1; plot!( msol.t, reduce(hcat, msol.u)'[:,i], color=[1 1] , alpha=0.75, lw=5 ) 
  i = 6; plot!( msol.t, reduce(hcat, msol.u)'[:,i], color=[1 1] , alpha=0.75, lw=5 ) 

  plot!( msol ; legend=false ) 

  b=[0.45, 0.4]
  
  external_forcing = rand(length(survey_time),6) # random
  efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
  hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
 
  p = ( b, K, d, dh, v, tau, hsa)   

  prob2 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=lags )
  # prob2 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )

  plot(0)
  msol2 =  solve( prob2,  solver, saveat=dt )
  i = 1; plot!( msol2.t, reduce(hcat, msol2.u)'[:,i], color=[3 3] , alpha=0.75, lw=5 ) 
  i = 6; plot!( msol2.t, reduce(hcat, msol2.u)'[:,i], color=[1 1] , alpha=0.75, lw=5 ) 
 
  plot!( msol2, label="dde, with random hsa" )

  plot!(; legend=true, xlim=(1990,2021) )


