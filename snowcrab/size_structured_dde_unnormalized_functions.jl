using Turing

function size_structured_dde!( du, u, h, p, t)
  # here u, du are actual numbers .. not normalized by K due to use of callbacks

  (b, K, d, d2, v, tau, hsa)  = p

  u1 = h(p, t-1.0)[2:5]    # no in previous years
  f8 = h(p, t-8.0)[6]  # no mature fem 8  yrs ago
  vh = hsa(t, 1:6)

  # this break down seems to speed it up a bit ... not sure why
  br =  f8 .* b
  tr =  v .* u1
  # dr =  d .* ( max.(u,0.0) ./ vh) .^ (2.0)
  uv = u ./ (K .* vh)
#  dr =  d .* u .+  d2 .* u .* uv
#  dr =  d .* u .+  d2 .* u .* ( uv .^ 2.0 )
dr =  d .* u .* uv .+  d2 .* u .* ( uv .^ 2.0 )
# dr =  d .* u .* uv .* (1.0 .+ uv )
 
  du[1] = tr[1]            - dr[1]       # note:
  du[2] = tr[2]   - tr[1]  - dr[2]
  du[3] = tr[3]   - tr[2]  - dr[3]
  du[4] = tr[4]   - tr[3]  - dr[4]
  du[5] = br[1]   - tr[4]  - dr[5]
  du[6] = br[2]            - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers

end

 

function dde_parameters()
    # these are dummy initial values .. just to get things started
    b=[2.1031608398813523, 2.957467881391256]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2 = d
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    params = ( b, K, d, d2, v, tau, hsa)
    return params
end
 

function affect_fishing!(integrator)
  i = findall(t -> t == integrator.t, fish_time)[1]
  integrator.u[1] -=  removed[ i[1] ]
  # integrator.u[1] -=  removed[ i ] / integrator.p[2][1]  # vintegrator.p[2] == K ; so divide by K[1]  
end



@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM,
    solver=MethodOfSteps(Tsit5()), dt = 0.01,  ::Type{T} = Float64) where T

    # biomass process model:
    K ~ filldist( TruncatedNormal( kmu, kmu*0.25, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all
 
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.2, 5.0), nS )
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -1.0, 1.0), nS )

    model_sd ~ truncated( Cauchy( 1.0, 1.0), 0.0, 1.0 )

    # birth rate from F(y - 8 to 10)  and for males m5 and females
    b ~ filldist( TruncatedNormal(10.0, 1.0, 0.01, 100.0), 2 )

    # background mortality
    d ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )
    d2 ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

    # transition rates
    v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 10.0), 4 )

    # initial conditions
    u0 ~ filldist( TruncatedNormal(0.6, 0.1, 0.01, 0.99), nS )

    pm = ( b, K, d, d2, v, tau, hsa )
    # @show pm

    # m = TArray{T}( length(Si), nS )

    # process model
    msol = solve(
        remake( prob; u0=u0 .* K, h=h, tspan=tspan, p=pm ),
        solver,
        callback=cb,
        # maxiters=1e6,
        # isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) , # force solutions to remain within 95% quantile bounds
        isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),
        saveat=dt
    )

    # @show msol.retcode
    if msol.retcode != :Success
      # fill this as without causes Turing confusion
      # for i in 1:nSI
      #     for k in 1:nS
      #         m[i,k] ~ TruncatedNormal( smallnumber, smallnumber , 0.0, 1.0 )  # dummy values as m has not been sampled yet .. without this it will crash
      #     end
      # end
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood of the data
    for i in 1:nSI
        ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
        for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] / K[k] * q[k] + qc[k], model_sd )  # observation and process error combined .. expand truncation for edge cases
        end

    end

end




function fishery_model_test( test=("basic", "random_external_forcing", "fishing", "nofishing") )

  ## test DifferentialEquations DDE model 
  gr()
  theme(:default)

  if any( occursin.( r"basic", test )  )

    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    b=[200, 100.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    d=[0.1, 0.2, 0.3, 0.3, 0.3, 0.3];
    d2 =d
    v=[0.9, 0.8, 0.8, 0.8];
    u0 = [ 0.65, 0.6, 0.52, 0.62, 0.58, 0.32 ] .* K ;
    tau=[1.0] 

    params = ( b, K, d, d2, v,  tau, hsa)
    prob1 = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]
    out = msol1 =  solve( prob1,  solver, callback=cb, saveat=dt )
    pl = plot()
    pl = plot( pl, msol1, ; legend=false, xlim=(1999,2021), title="basic" )

  elseif  any( occursin.( r"random_external_forcing", test )  )

    # testing DDE version -- initial values for size_structured_dde! 
    # basic run
    ks = 1000
    u0 = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ]  .* ks;
    b=[ 10, 10 ]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    d2 =d
    v=[0.9, 0.9, 0.9, 0.9];  
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
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks .*0.5
    tau=[1.0] 
    
    solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
    prob2 = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=[tau] )
    out = msol2 =  solve( prob2,  solver, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol2 ; title="random_external_forcing" )
  
  elseif  any( occursin.( r"nofishing", test ) )

    ks = 1.0e9
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
        
    # K=[1.1309624365561053e7, 6.65896205845559e6, 1.1175886847805813e7, 8.096057237318473e6, 4.345179415118201e6, 7.880503154236994e6]
    u0 = K .*0.75;
     b=[2.1031608398813523, 2.957467881391256]
     d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
     d2 =d
     v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks
  
    tau=[1.0] 
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    external_forcing = ones(length(survey_time),6)  # turn off external forcing  (wrt viable habitat)
    efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
    p = ( b, K, d, d2, v, tau, hsa )   
  
    prob3 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=[tau] )
    out = msol3 =  solve( prob3,  solver, saveat=dt  ) #, isoutofdomain=(y,p,t)->any(x->(x<0)|(x>1), y) )
    pl = plot()

    pl = plot!( pl, msol3, title="dde, with random hsa, with fishing") 

    i = 1; # fb
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_1" ) 

    i = 6; # fem
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_6" ) 
  
   
  elseif  any( occursin.( r"fishing", test )  | occursin( r"nofishing", test ) )

    ks = 1.0e9
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = K .*0.75;
    b=[2.1031608398813523, 2.957467881391256]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2 =d
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks

    tau=[1.0] 
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    # add random viable habitat forcing
    external_forcing = rand(length(survey_time),6) # random
    efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
   
    p = ( b, K, d, d2, v, tau, hsa)   
  
    prob4 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=[tau] )
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
    z > nZ && break
    b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
    d2 = [ res[j, Symbol("d2[$k]"), l] for k in 1:nS]

    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
    qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]
    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

    pm = ( b, K, d, d2, v, tau, hsa )

    prb = remake( prob; u0=u0.*K , h=h, tspan=tspan, p=pm )

    msol = solve( prb, solver, callback=cb, saveat=dt )
    msol2 = solve( prb, solver, saveat=dt ) # no call backs

    for i in 1:nM
        ii = findall(x->x==prediction_time[i], msol.t)[1]
        jj = findall(x->x==prediction_time[i], msol2.t)[1]
        sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t[ii])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
        sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t[jj]) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
        md[i,:,z,1] = msol.u[ii]  ./ K # with fishing
        md[i,:,z,2] = msol2.u[jj]  ./ K # no fishing
        mn[i,:,z,1] = msol.u[ii]     # with fishing
        mn[i,:,z,2] = msol2.u[jj]   # no fishing
        mb[i,z,1] = mn[i,1,z,1]  .* sf
        mb[i,z,2] = mn[i,1,z,2]  .* sf2

    end

  end
  end

  # extract sims (with fishing)
  g = mb[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
  # size(g)

  # plot biomass
  gr()
  pl = plot()
  pl = plot!(pl, prediction_time, g[:,sample(1:nZ, nI)];  alpha=0.05, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  
  ub = quantile(vec(g), 0.95)
  pl = plot!(pl; ylim=(0, ub ) )

  # back transform S to normal scale .. do sims too (TODO)
  k = 1
  yhat = ( S[:,k] .- mean(res[:,Symbol("qc[$k]"),:]  ) ) ./ mean(res[:,Symbol("q[$k]"),:]) .* mean(res[:,Symbol("K[$k]"),:]  )
  # yhat = abundance_from_index( S[:,k], res, k ) 


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


function fishery_model_predictions_trace( res; n_sample=10, plot_k=1, alpha=0.01  )

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
        d2 = [ res[j, Symbol("d2[$k]"), l] for k in 1:nS]

        v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]

        q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
        qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]

        u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

        pm = ( b, K, d, d2, v, tau, hsa )

        prb = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

        if plot_k==1
            # do fishing and nonfishing

            msol = solve( prb, solver, callback=cb, saveat=dt )
            msol2 = solve( prb, solver, saveat=dt ) # no call backs

            sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t) ./ 1000.0 ./ 1000.0 :  scale_factor
            sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt

            yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k])   .* sf2
            yval = vec( reduce(hcat, msol.u)'[:,plot_k] )   .* sf

            pl = plot!( pl, msol.t, yval, alpha=alpha, lw=1, color=:orange )
            pl = plot!( pl, msol2.t, yval2, alpha=alpha*2.0, lw=1, color=:lime )

            push!(out, yval)
            push!(out2, yval2)

        else
            msol2 = solve( prb, solver, saveat=dt ) # no call backs
            yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k])  
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

