
using Turing

# m = abundance (prefishery)
# survey timing: spring before 2004 and fall afterwards 
# fishery operates in winter for 4x and spring and summer for N and Sens
#   4x:  (m[i-1] - removals[i-1]) = abundance (post fishery) upon which dynamics is applied, to give m[i], the abundance (prefishery), so S[i] ~ m[i] - rem[i]
#   n and sens: post 2004 .. same as above, , so S[i] ~ m[i] - rem[i] (surveys post fishery)
#               pre 2004 ... same as above but.. S[i] ~ m[i]  (no removals due to survey being prefishery)

@model function logistic_discrete_turing( S, kmu, nT, nM, removed, ::Type{T} = Float64) where T 
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  

  K ~ TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.5 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.5 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
  qc ~ TruncatedNormal( SminFraction, 0.1, -1.0, 1.0) 

  m = TArray{T}( nM )
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
  end

  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
 
  # likelihood
  # observation model: Y = q X + qc ; X = (Y - qc) / q
  for i in iok
    S[i] ~ Normal( q * (m[i] - removed[i]/K) + qc, bosd )  ;
  end 

end




@model function logistic_discrete_turing_basic( S, kmu, nT, nM, removed, ::Type{T} = Float64) where T
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.5, 1.5)    

  m = TArray{T}( nM )
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
  end
  
  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ; # predict with no removals
  end
    
  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   
  # likelihood
  # observation model: Y = q X  ; X = (Y ) / q
   for i in iok
    S[i] ~ Normal( q * (m[i] - removed[i]/K), bosd )  ;
  end 

end
  

@model function logistic_discrete_turing_historical( S, kmu, nT, nM, removed, ::Type{T} = Float64) where T
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~  TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0 )   
  r ~  TruncatedNormal( 1.0, 0.1, 0.25, 2.0)   # (mu, sd)

  bpsd ~  truncated( Cauchy( 0, 0.1), 1.0e-9, 0.5 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated( Cauchy( 0, 0.1), 1.0e-9, 0.5 )    ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  1.0e-9, 10.0 )    

  # m's are "total avaialble for fishery"
  m = TArray{T}( nM )
  m[1] ~ truncated( Beta( 8, 2) )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.25)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.25)  ; # predict with no removals
  end
 
  if any( x -> x < 0.0, m)  # permit overshoot
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood
  # observation model: Y = q X  ; X = (Y ) / q
  s = S ./ K
  for i in iok
    s[i] ~ Normal( q * ( m[i] - removed[i]/K ) , bosd )  ; # fall survey
  end

end
  


@model function logistic_discrete_map_turing( S, kmu, nT, nM, removed,  ::Type{T} = Float64) where T
  # biomass process model: n(t+1) = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

    K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
    r ~  TruncatedNormal( 1.0, 0.1, 0.5, 3.0)   # (mu, sd)

    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
    qc ~ TruncatedNormal( 0.0, 0.1, -1.0, 1.0) 

    m = TArray{T}( nM )# fished (postfishery) abundance
    m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

    for i in 2:nT
      m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
    end

    for i in (nT+1):nM
      m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
    end

    if any( x -> x < 0.0, m)  # permit overshoot
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    # likelihood
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for i in iok
      S[i] ~ Normal( q * (m[i] - removed[i-1]/K ) + qc, bosd )  ;
    end
  
end


@model function logistic_discrete_turing_north_south( S, kmu, nT, nM, removed, ty = 6, ::Type{T} = Float64) where T 
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 
  # NENS and SENS

  K ~ TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.5 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.5 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
  qc ~ TruncatedNormal( SminFraction, 0.1, -1.0, 1.0) 

  m = TArray{T}( nM )
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
  end

  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
 
  # likelihood
  # observation model: Y = q X + qc ; X = (Y - qc) / q
   
  for i in iok
    if i < ty
      S[i] ~ Normal( q * ( m[i] ) + qc, bosd )  ;  # spring survey
    elseif i == ty
      S[i] ~ Normal( q * ( m[i] - (removed[i-1] + removed[i]) / (2.0*K ) )+ qc, bosd )  ;  # transition year 
    else
      S[i] ~ Normal( q * ( m[i] - removed[i]/K )+ qc, bosd )  ; # fall survey
    end
  end

end



@model function logistic_discrete_turing_basic_north_south( S, kmu, nT, nM, removed, ty = 6, ::Type{T} = Float64) where T
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.5, 1.5)    

  m = TArray{T}( nM )
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
  end
  
  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ; # predict with no removals
  end
  
  
  if any( x -> x < 0.0 || x >1.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   
  # likelihood
  # observation model: Y = q X  ; X = (Y ) / q
   
  for i in iok
    if i < ty
      S[i] ~ Normal( q * ( m[i] ) , bosd )  ;  # spring survey
    elseif i == ty
      S[i] ~ Normal( q * ( m[i] - (removed[i-1] + removed[i]) / (2.0*K ) ) , bosd )  ;  # transition year 
    else
      S[i] ~ Normal( q * ( m[i] - removed[i]/K ) , bosd )  ; # fall survey
    end
  end

end
  
 
@model function logistic_discrete_turing_historical_north_south( S, kmu, nT, nM, removed, ty = 6, ::Type{T} = Float64) where T
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~  TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0 )   
  r ~  TruncatedNormal( 1.0, 0.1, 0.25, 2.0)   # (mu, sd)

  bpsd ~  truncated( Cauchy( 0, 0.1), 1.0e-9, 0.5 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  truncated( Cauchy( 0, 0.1), 1.0e-9, 0.5 )     ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  1.0e-9, 10.0 )    

  # m = total available biomass
  m = TArray{T}( nM )
  m[1] ~  truncated( Beta( 8, 2) )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.25)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.25)  ; # predict with no removals
  end

  if any( x -> x < 0.0 || x >1.25, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   
  # likelihood
  # observation model: Y = q X  ; X = (Y ) / q

  # yrs = 1999:2021  # <<<<<<<<-- change
  # spring to fall survey: transition year = 2004
  # spring = 1:5
  # fall = 6:last
  # in cfa4x fishery always after survey
  # m's are "postfishery"
  
  s = S ./ K
  
  for i in iok

    if  i < ty
      s[i] ~ Normal( q * ( m[i] ), bosd )  ;  # spring survey
    elseif i == ty
      s[i] ~ Normal( q * ( m[i] - (removed[i-1] + removed[i]) / (K*2.0 ) ), bosd )  ;  # transition year 
    else
      s[i] ~ Normal( q * ( m[i] - removed[i] /K), bosd )  ; # fall survey
    end
  end
    

end
  


  
function fishery_model_test( test=("basic" ) )

  ## test model by sampling from random priors 
  gr()
  theme(:default)
  pl = plot()

  if any( occursin.( r"basic", test )  )

    res = sample( fmod, Prior(), 100, nwarmup = 100, nchains =1 )
 
    for l in 1:size(res)[3]
      for i in 1:length(res)  
          w = zeros(nM)
          for j in 1:nM
              w[j] = res[i, Symbol("K"),l] * res[i, Symbol("m[$j]"),l] 
          end
          pl = plot!(pl, prediction_time, w;  alpha=0.1, color=:orange)
      end
    end
    pl = plot!(pl; legend=false, title="basic prior check" )
 
  end
 
  return (res, pl) 

end



function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=100 )
  # n_sample = num samples to plot
  nchains = size(res)[3]
  nsims = size(res)[1]

  nZ = nchains*nsims
  nI = Int( min( nZ , n_sample ) )
 
  mb = md = zeros(nM, nZ)  # biomass normalized
  
  z = 0
  for j in 1:nsims  # nsims
  for l in 1:nchains #nchains
    z += 1
    for i in 1:nM
        md[i,z] = res[j, Symbol("m[$i]"), l]
        mb[i,z] = md[i,z]  * res[j, Symbol("K"), l]
    end
  end
  end

  # plot biomass
  gr()
  pl = plot()
  pl = plot!(pl, prediction_time, mb[:,sample(1:nZ, nI)];  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(mb, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(mb)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  # S[i,k] ~ TruncatedNormal( msol.u[ii][k] * q[k] + qc[k], bpsd, 0.0, 1.0)
  if model_variation=="logistic_discrete_basic"  
    yhat = S ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  )
  elseif model_variation=="logistic_discrete_historical"  
    yhat = ( S ) ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  )
  elseif model_variation=="logistic_discrete"  
    yhat = ( S .- mean(res[:,Symbol("qc"),:] ) ) ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  )
  elseif model_variation=="logistic_discrete_map"  
    yhat = ( S .- mean(res[:,Symbol("qc"),:] ) ) ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  )
  end

  pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )

  return (md, nothing, mb, pl)

end


# -----------


function fishery_model_predictions_timeseries(nothing; prediction_time=prediction_time, plot_k=1)
  print("Nothing to do in model with a single state variable")
  return (nothing, nothing)
end


# ----------

function fishery_model_predictions_trace( res; n_sample=10, plot_k=1, alpha=0.01, plot_only_fishing=true )
  print("Nothing to do in model with a single state variable")
  return (nothing, nothing, nothing)
end


# ----------

function logistic_discrete_reference_points(r, K)
  expK = exp.(K) 
  msy   = r .* expK ./ 4.0 ; # maximum height of of the latent productivity (yield)
  bmsy  = expK ./ 2.0 ; # biomass at MSY
  fmsy  = 2.0 .* msy ./ expK ; # fishing mortality at MSY
  return (msy, bmsy, fmsy)
end

      # plot!(prediction_time, u, lw=2, color=:orangered )
      # scatter!(prediction_time, u, markersize=4, color=:goldenrod1 )

      # plot!(survey_time, yhat, color=:purple2, lw=2 )
      # scatter!(survey_time, yhat, markersize=4, color=:purple4)
 
function fishery_model_harvest_control_rule(res, yrs; FM=FM, fb=fb, n_sample=500 )

  fmsy = nothing

  pl = plot()
 
  r = vec( Array(res[:, Symbol("r"), :]) )
  K = vec( Array(res[:, Symbol("K"), :]) ) 
  (msy, bmsy, fmsy) = logistic_discrete_reference_points(r, K)
  pl = hline!(pl, sample(fmsy, n_sample); alpha=0.01, color=:lightgray )
  pl = hline!(pl, [mean(fmsy)];  alpha=0.6, color=:darkgray, lw=5 )
  pl = hline!(pl, [quantile(fmsy, 0.975)];  alpha=0.5, color=:gray, lw=2, line=:dash )
  pl = hline!(pl, [quantile(fmsy, 0.025)];  alpha=0.5, color=:gray, lw=2, line=:dash )


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
  
  Fkt = removed
 
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

