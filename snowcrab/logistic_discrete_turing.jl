
using Turing 


if model_variation == "Logistic_q"

  @model function logistic_discrete_turing( S, kmu, nT, nM, removed, solver=MethodOfSteps(Tsit5()) )
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors 
  
    K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
    r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)
  
    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
    q ~ TruncatedNormal(  1.0, 0.1,  0.5, 1.5)    
  
    m =  tzeros( nM ) # fished (postfishery) abundance
    m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event
  
    for i in 2:nT
      m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
    end
    
    for i in (nT+1):nM
      m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ; # predict with no removals
    end
    
    
    if any( x -> x < smallnumber, m)
      Turing.@addlogprob! -Inf
      return nothing
    end
    
    iok = findall( !ismissing, S )
  
    # likelihood
    # observation model: Y = q X  ; X = (Y ) / q
    for i in iok
      S[i] ~ TruncatedNormal( q * m[i], bosd, 0.0, 1.0 )  ;
    end
  end
  
elseif model_variation == "Logistic_q_qc"

  @model function logistic_discrete_turing( S, kmu, nT, nM, removed, iok=findall( !ismissing, S ), 
    solver=MethodOfSteps(Tsit5()) ) 
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors 

    K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
    r ~  TruncatedNormal( 1.0, 0.1, 0.5, 1.5)   # (mu, sd)

    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
    qc ~ TruncatedNormal( 0.0, 0.1, -0.5, 0.5) 

    m =  tzeros( nM ) # fished (postfishery) abundance
    m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event

    for i in 2:nT
      m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
    end

    for i in (nT+1):nM
      m[i] ~ TruncatedNormal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
    end

    if any( x -> x < smallnumber, m)
      Turing.@addlogprob! -Inf
      return nothing
    end

    #  check positivity of back transform
    yhat = S[iok] .*  q  .-  qc   
    if any( x -> x < smallnumber, yhat )
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for i in iok
      S[i] ~ TruncatedNormal( q * m[i] + qc, bosd, 0.0, 1.0 )  ;
    end
  end

elseif model_variation == "Logistic_Map"
  
  @model function logistic_discrete_turing( S, kmu, nT, nM, removed, solver=MethodOfSteps(Tsit5()) )
    # biomass process model: n(t+1) = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors 
  
      K ~ TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0)  
      r ~  TruncatedNormal( 1.0, 0.1, 0.5, 3.0)   # (mu, sd)
  
      bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
      bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
      q ~ TruncatedNormal(  1.0, 0.1,  0.01, 10.0)    
      qc ~ TruncatedNormal( 0.0, 0.1, -0.5, 0.5) 
  
      m =  tzeros( nM ) # fished (postfishery) abundance
      m[1] ~  TruncatedNormal( 0.9, 0.2, 0.1, 1.0 )  ; # starting b prior to first catch event
  
      for i in 2:nT
        m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 0.0, 1.0)  ;
      end
  
      for i in (nT+1):nM
        m[i] ~ TruncatedNormal(  r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 0.0, 1.0)  ;  # predict with no removals
      end
  
      if any( x -> x < smallnumber, m)
        Turing.@addlogprob! -Inf
        return nothing
      end
  
      #  check positivity of back transform
  
      yhat = S[iok] .*  q  .-  qc   
      if any( x -> x < smallnumber, yhat )
        Turing.@addlogprob! -Inf
        return nothing
      end
  
      # likelihood
      # observation model: Y = q X + qc ; X = (Y - qc) / q
      for i in iok
        S[i] ~ TruncatedNormal( q * m[i] + qc, bosd, 0.0, 1.0 )  ;
      end
  end
  
  
else

  print( "model_variation is not correct?" )

end




  

 
