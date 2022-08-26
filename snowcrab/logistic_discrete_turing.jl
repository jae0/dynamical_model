
using Turing 

@model function logistic_discrete_turing( S, kmu, nT, nP, removed, solver=MethodOfSteps(Tsit5()), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors 

    K ~ TruncatedNormal( kmu, kmu*0.1, kmu/2.5, kmu*2.5)  
    r ~  TruncatedNormal( 1.0, 0.1, 0.25, 1.75)   # (mu, sd)

    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

    q ~ TruncatedNormal(  1.0, 0.1,  0.75, 1.5)    
    qc ~ TruncatedNormal( 0.0, 0.1, -1.0, 1.0) 
 
    nM = nT + nP

    m =  Vector{T}(undef, nM ) # fished (postfishery) abundance
    m[1] ~  TruncatedNormal( 0.9, 0.2, 0.5, 1.25 )  ; # starting b prior to first catch event

    for i in 2:nT
      m[i] ~ TruncatedNormal( r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 1.0e-9, 1.25)  ;
    end

    for i in (nT+1):nM
      m[i] ~ TruncatedNormal( r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 1.0e-9, 1.25)  ;
    end

    if any( x -> x < eps, m)
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood
    for i in 1:nT
      S[i] ~ TruncatedNormal( (m[i] + qc)*q, bosd, 1.0e-9, 1.25 )  ;
    end
  end



@model function logistic_discrete_turing_basic( S, kmu, nT, nP, removed, solver=MethodOfSteps(Tsit5()), ::Type{T}=Float64 ) where {T}  
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 

  K ~ TruncatedNormal( kmu, kmu*0.1, kmu/2.5, kmu*2.5)  
  r ~  TruncatedNormal( 1.0, 0.1, 0.25, 1.75)   # (mu, sd)

  bpsd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  bosd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)

  q ~ TruncatedNormal(  1.0, 0.1,  0.75, 1.5)    
 
  nM = nT + nP
  m =  Vector{T}(undef, nM) # fished (postfishery) abundance
  m[1] ~  TruncatedNormal( 0.9, 0.2, 0.5, 1.25 )  ; # starting b prior to first catch event

  for i in 2:nT
    m[i] ~ TruncatedNormal( r * m[i-1] * ( 1.0 - m[i-1] ) - removed[i-1]/K, bpsd, 1.0e-9, 1.25)  ;
  end

  for i in (nT+1):nM
    m[i] ~ TruncatedNormal( r * m[i-1] * ( 1.0 - m[i-1] ), bpsd, 1.0e-9, 1.25)  ;
  end
  
  if any( x -> x < eps, m)
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood
  for i in 1:nT
    S[i] ~ TruncatedNormal( m[i]*q, bosd, 1.0e-9, 1.25 )  ;
  end
end
