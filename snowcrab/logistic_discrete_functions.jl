using Turing

@model function logistic_discrete_turing( S, kmu, nT, nM, removed, ::Type{T} = Float64) where T 
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
  qc ~ TruncatedNormal( 0.0, 0.1, -0.5, 0.5) 
 
  m = TArray{T}( nM )
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
  end

  if any( x -> x < 0.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   

  #  check positivity of back transform
  # yhat = S[iok] .*  q  .-  qc   
  # if any( x -> x < 0.0, yhat )
  #   Turing.@addlogprob! -Inf
  #   return nothing
  # end

  # likelihood
  # observation model: Y = q X + qc ; X = (Y - qc) / q
  for i in iok
    S[i] ~ TruncatedNormal( q * m[i] + qc, bosd, 0.0, 1.0 )  ;
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
  
  
  if any( x -> x < 0.0, m)
    Turing.@addlogprob! -Inf
    return nothing
  end
   
  # likelihood
  # observation model: Y = q X  ; X = (Y ) / q
  for i in iok
    S[i] ~ TruncatedNormal( q * m[i], bosd, 0.0, 1.0 )  ;
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
    qc ~ TruncatedNormal( 0.0, 0.1, -0.5, 0.5) 

    m = TArray{T}( nM )# fished (postfishery) abundance
    m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

    for i in 2:nT
      m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
    end

    for i in (nT+1):nM
      m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
    end

    if any( x -> x < 0.0, m)
      Turing.@addlogprob! -Inf
      return nothing
    end

    # #  check positivity of back transform
    # yhat = S[iok] .*  q  .-  qc   
    # if any( x -> x < 0.0, yhat )
    #   Turing.@addlogprob! -Inf
    #   return nothing
    # end

    # likelihood
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for i in iok
      S[i] ~ TruncatedNormal( q * m[i] + qc, bosd, 0.0, 1.0 )  ;
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



function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=1e10 )
 
  nchains = size(res)[3]
  nsims = size(res)[1]

  nI = Int( min( nchains*nsims , n_sample ) )

  mb = md = zeros(nM, nI)  # biomass normalized
  
  z = 0

  for j in 1:nsims  # nsims
  for l in 1:nchains #nchains
    z += 1
    z > nI && break
    for i in 1:nM
        md[i,z] = res[j, Symbol("m[$i]"), l]
        mb[i,z] = md[i,z]  * res[j, Symbol("K"), l]
    end
  end
  end

    # plot biomass
  gr()
  pl = plot()
  pl = plot!(pl, prediction_time, mb;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(mb, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(mb)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  # S[i,k] ~ TruncatedNormal( msol.u[ii][k] * q[k] + qc[k], bpsd, 0.0, 1.0)
  if model_variation=="logistic_discrete_basic"  
    yhat = S ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  )
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
 