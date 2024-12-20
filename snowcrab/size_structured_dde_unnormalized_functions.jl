using Turing

function size_structured_dde!( du, u, h, p, t)
  # here u, du are actual numbers .. not normalized by K due to use of callbacks

  (b, K, d, d2, v, tau, hsa)  = p

  # u1 = h(p, t-1.0)[2:5]    # no in previous years
  # f8 = h(p, t-8.0)[6]      # no mature fem 8  yrs ago
  uu = u ./ K 
  uh = uu ./ hsa(t, 1:6)          # viable habitat surface area 

  # this break down seems to speed it up a bit ... not sure why
  br =  b .* h(p, t-8.0)[6]  
  tr =  v .* h(p, t-1.0)[2:5]

  # d: background mortality + d2: excess mortality due to habitat fluctuations 
  dr =  d .* u .+  (d2 .*  uh ) .* uh  #  death rate

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

 

@model function size_structured_dde_turing( S, kmu, tspan, prob,  nS,  
    solver=MethodOfSteps(Tsit5()), dt = 0.01, ::Type{T} = Float64) where T

    # biomass process model:
    # K ~ filldist( TruncatedNormal( kmu, kmu*0.2, kmu/1000.0, kmu*1000.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all
    K ~ filldist( truncated( LogNormal(log(kmu), 0.25), kmu/1000, kmu*1000), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.2, 5.0), nS )
    qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 

    model_sd ~ filldist( Gamma( 2.0, 0.1 ), nS)

    # birth rate from F(y - 8 to 10)  and for males m5 and females
    # b ~ filldist( TruncatedNormal(5.0, 0.1, 0.01, 100.0), 2 )
    b ~ filldist( truncated(Chisq( 3 ), 0.1, 10),  2 )   # centered on 1; plot(x->pdf(Chisq(7), x), xlim=(0,10)) # mode of 5

    # background mortality
    d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 1.0), nS )

    # excess mortality due to habitat (second order)
    d2 ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 1.0), nS )

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.01, 1.0), 4 )

    # initial conditions
    u0 ~ filldist( TruncatedNormal(0.5, 0.3, 0.01, 1.0), nS )

    pm = ( b, K, d, d2, v, tau, hsa )
    # @show pm

    # m = TArray{T}( length(Si), nS )

    # process model
    msol = solve(
        remake( prob; u0=u0 .* K, h=h, tspan=tspan, p=pm ),
        solver,
        callback=cb,
        # maxiters=1e5,
        abstol = 1e-9, 
        # reltol = 1e-9
        isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),
        saveat=dt
    )

    # @show msol.retcode
    if !SciMLBase.successful_retcode(msol)  
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood of the data
    for i in 1:nSI
        ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
        for k in 1:nS
          S[Si[i],k] ~ Normal( ( msol.u[ii][k] ./ K[k] .* q[k] + qc[k]) .* K[k], model_sd[k] )  # observation and process error combined .. expand truncation for edge cases
        end

    end

end




function fishery_model_test( test=("basic", "random_external_forcing", "fishing", "nofishing") )

  ## test DifferentialEquations DDE model 
  gr()
  theme(:default)

  if any( occursin.( r"basic", test )  )

    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    b=[3, 3]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    d=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
    d2 = [0.1, 0.2, 0.2, 0.2, 0.2, 0.2]
    v=[0.6, 0.7, 0.6, 0.7];
    u0 = [ 0.65, 0.6, 0.52, 0.62, 0.58, 0.32 ] .* K ;
    tau=[1.0] 
    
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* kmu

    survey_time = 1999:2021
    tspan = (1990.0, 2021.0)
    dt = 0.001
    
    external_forcing = ones(length(survey_time),6)  # turns it off
    # external_forcing = rand(length(survey_time),6) # random
    
    efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc, 1999:2021, 1:6 )
    solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
 
    params = ( b, K, d, d2, v, tau, hsa)
    prob1 = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]
    out = msol1 =  solve( prob1,  solver, callback=cb, saveat=dt )
    pl = plot( msol1, ; legend=false, xlim=(1999,2021), title="basic" )

  elseif  any( occursin.( r"random_external_forcing", test )  )

    # testing DDE version -- initial values for size_structured_dde! 
    # basic run
    ks = 1000
    u0 = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ]  .* ks;
    b=[ 1, 1 ]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    d2=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
    d2=[0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
    v=[0.9, 0.9, 0.9, 0.9];  
    tau=[1.0] 
  
    survey_time = 1999:2021
    
    external_forcing = ones(length(survey_time),6)  # turns it off
    # external_forcing = rand(length(survey_time),6) # random
    
    efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc, 1999:2021, 1:6 )
  
    p = ( b, K, d, d2, v, tau, hsa )   
    tspan = (1990.0, 2010.0)
    nS = length(u0)  # n components
    
    # history function for time before t0:: 0.5 default
    # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks .*0.5
    tau=[1.0] 
    
    solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
    prob2 = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=tau )
    out = msol2 =  solve( prob2,  solver, saveat=dt )
    pl = plot( msol2 ; title="random_external_forcing" )
  
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

  

function fishery_model_predictions_old( res; prediction_time=prediction_time, n_sample=100 )

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

    msol = solve( prb, solver, callback=cb, saveat=dt  )
    msol2 = solve( prb, solver, saveat=dt   ) # no call backs
    
    if (SciMLBase.successful_retcode(msol)) && (SciMLBase.successful_retcode(msol2)) 

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


function fishery_model_predictions_trace_old( res; n_sample=10, plot_k=1, alpha=0.01  )

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

            msol = solve( prb, solver, callback=cb, saveat=dt  )
            msol2 = solve( prb, solver, saveat=dt  ) # no call backs

            if (!SciMLBase.successful_retcode(msol)  ) && (!SciMLBase.successful_retcode(msol2)  ) 

              sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t) ./ 1000.0 ./ 1000.0 :  scale_factor
              sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt

              yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k])   .* sf2
              yval = vec( reduce(hcat, msol.u)'[:,plot_k] )   .* sf

              pl = plot!( pl, msol.t, yval, alpha=alpha, lw=1, color=:orange )
              pl = plot!( pl, msol2.t, yval2, alpha=alpha*2.0, lw=1, color=:lime )

              push!(out, yval)
              push!(out2, yval2)
            end
        else
          if (!SciMLBase.successful_retcode(msol)  ) && (!SciMLBase.successful_retcode(msol2)  ) 

            msol2 = solve( prb, solver, saveat=dt ) # no call backs
            yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k])  
            pl = plot!( pl, msol2.t, yval2, alpha=alpha, lw=1, color=:lime )

            push!(out2, yval2)
          end
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
  colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]

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



# -------------------


function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=100 )

  nchains = size(res)[3]
  nsims = size(res)[1]
  
  md = zeros(nM, nS, n_sample, 2)  # number normalized
  mn = zeros(nM, nS, n_sample, 2)  # numbers
  mb = mn[:,1,:,:]  # biomass of first class

  trace_time = Vector{Vector{Float64}}()

  out10 = Vector{Vector{Float64}}()
  out11 = Vector{Vector{Float64}}()
  out12 = Vector{Vector{Float64}}()
  
  out1 = Vector{Vector{Float64}}()
  out2 = Vector{Vector{Float64}}()
  out3 = Vector{Vector{Float64}}()
  out4 = Vector{Vector{Float64}}()
  out5 = Vector{Vector{Float64}}()
  out6 = Vector{Vector{Float64}}()

  ntries = 0
  z = 0

  while z <= n_sample 
    ntries += 1
    ntries > n_sample*10 && break
    z >= n_sample && break

    j = rand(1:nsims)  # nsims
    l = rand(1:nchains) #nchains
    b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
    d2=[ res[j, Symbol("d2[$k]"), l] for k in 1:nS]
    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

    pm = ( b, K, d, d2, v, tau, hsa )
    prb = remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm )
    msol1 = solve( prb, solver, callback=cb, saveat=dt, dt=dt  )

    z==0 && push!(trace_time, msol1.t)

    msol0 = solve( prb, solver, saveat=trace_time[1] ) # no call backs
    
    if !SciMLBase.successful_retcode(msol1)   && !SciMLBase.successful_retcode(msol0)  
        z += 1
        
        # annual
        for i in 1:nM
            ii = findall(x->x==prediction_time[i], trace_time[1]) 
            if length(ii) > 0  
              ii = ii[1]
              sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol1.t[ii])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
              md[i,:,z,1] = msol1.u[ii]  ./ K # with fishing
              md[i,:,z,2] = msol0.u[ii]  ./ K # no fishing
              mn[i,:,z,1] = msol1.u[ii]   # with fishing scaled to K
              mn[i,:,z,2] = msol0.u[ii]   # no fishing scaled to K
              mb[i,z,1] = mn[i,1,z,1]  .* sf  # biomass of state var 1 
              mb[i,z,2] = mn[i,1,z,2]  .* sf
            end
        end  # end for

        # traces for plotting, etc 
        sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(trace_time[1])  ./ 1000.0 ./ 1000.0 :  scale_factor
        b10 = vec( reduce(hcat, msol0.u)'[:,1]) .* sf
        b11 = vec( reduce(hcat, msol1.u)'[:,1]) .* sf
        b12 = ( b10 .- b11 ) ./ b10 

        push!(out10, b10)
        push!(out11, b11)
        push!(out12, b12)

        push!(out1, vec( reduce(hcat, msol1.u)'[:,1]) )
        push!(out2, vec( reduce(hcat, msol1.u)'[:,2]) )
        push!(out3, vec( reduce(hcat, msol1.u)'[:,3]) )
        push!(out4, vec( reduce(hcat, msol1.u)'[:,4]) )
        push!(out5, vec( reduce(hcat, msol1.u)'[:,5]) )
        push!(out6, vec( reduce(hcat, msol1.u)'[:,6]) )

      end # if
  end  # while
 
  if z < n_sample 
    @warn  "Insufficient number of solutions" 
  end

  trace = (out1, out2, out3, out4, out5, out6 )
  trace_bio = (out10, out11, out12)
  return (md, mn, mb, trace, trace_bio, trace_time[1] )

end



# -----------




function fishery_model_mortality_old( removed, fb; n_sample=100 )    
  
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


function fishery_model_plot(; toplot=("fishing", "nofishing", "survey"),
  res=res, bio=bio, num=num, trace=trace, trace_bio=trace_bio, FM=FM, 
  S=S, si=1, scale_factor=scale_factor, 
  prediction_time=prediction_time, survey_time=survey_time, trace_time=trace_time, yrs=yrs, 
  alphav=0.075, pl= plot(), time_range=(floor(minimum(survey_time))-1.0, ceil(maximum(survey_time))+1.0 )
)
 
  # extract sims (with fishing)
  # plot biomass
  if any(isequal.("fishing", toplot))  
    g = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("nofishing", toplot))  
    g = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:lime)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("footprint", toplot))  
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:lightslateblue)
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
    yhat = ( S[:,k] ./ Km .- qcm  ) ./ qm .* Km   # abundance_from_index S[:,1]  
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
    pl = plot!(pl, prediction_time, gk;  alpha=alphav, color=:lightslateblue)
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


  if any(isequal.("trace", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl, trace_time, trace_bio[1], alpha=alphav, lw=1, color=:lime )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_fishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2], alpha=alphav, lw=1, color=:orange )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end
 
  if any(isequal.("trace_nofishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[1], alpha=alphav, lw=1, color=:lime )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end


  if any(isequal.("trace_footprint", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3], alpha=alphav, lw=1, color=:lime )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end

  if any(isequal.("fishing_mortality", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, survey_time, FM ;  alpha=0.02, color=:lightslateblue)
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
    pl = scatter!(pl, FM, g;  alpha=alphav, color=:lightslateblue)
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
  
    pl = vline!(pl, K;  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K./4;  alpha=0.05, color=:darkred )
  
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
  
    pl = vline!(pl, K;  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K./4;  alpha=0.05, color=:darkred )
  
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


 



function plot_prior_posterior( vn, prior, posterior; bw=0.02 )
  pri =  vec(collect( prior[:,Symbol(vn),:] ))
  pos =  vec(collect( posterior[:,Symbol(vn),:] ))
  pl = plot(ylabel="Density", xlabel=vn ) 
  pl = density!(pl, pri,  fill=true, color = :slateblue, fillalpha=0.25, bandwidth = bw, lw=0, label="Prior")
  pl = density!(pl, pos,  fill=true, color = :purple, fillalpha=0.5, bandwidth = bw, lw=0, label="Posterior")
  return(pl)
end

 
