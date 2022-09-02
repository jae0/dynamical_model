
using Turing, DifferentialEquations

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nP, solver=MethodOfSteps(Tsit5()) )   
    
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/5.0, kmu*5.0), nS )  
  
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -1.0, 1.0), nS )  
  
       
    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.001, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.1, 10.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.5, 0.1, 0.25, 0.9), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.5, 1.0), 4 ) 


    # initial conditions
    nM = nT + nP
    m = tzeros(nM, nS)
    for k in 1:nS 
      m[1,k] ~  TruncatedNormal( 0.9, 0.2, 0.5, 1.25 )  ; # starting b prior to first catch event permit higher than 1.0
    end 
    
    u0 = [ m[1,1], m[1,2], m[1,3], m[1,4], m[1,5], m[1,6]  ] .* K  # taking array slice keeps dual info .. but solve just wants floats .. reconstruct vector
 
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
     
    msol = solve( prob, 
      solver, 
      callback=cb, 
      saveat=dt, 
      isoutofdomain=(u,p,t)->any(x -> x < eps, u)
    ) 
     
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    for i in 2:nM
      j = findall(t -> t==prediction_time[i], msol.t)
      if length(j) > 0
        usk = max.( msol.u[j[1]], 1.0) ./ K 
        for k in 1:nS
          m[i,k] ~ TruncatedNormal( usk[k], bpsd, 1.0e-9, 1.25)  ; 
        end
      end
    end
  
    # in case no viable solutions predicetd at required time steps
    if any( x -> x < eps, m)
      Turing.@addlogprob! -Inf
      return nothing
    end

    # check positivity of back transform
    # column-wise arithmetic
    if any( x -> x < eps, (S .+ qc') .*  q' .*  K' )
      Turing.@addlogprob! -Inf
      return nothing
    end

    # observation model
    for k in 1:nS
      for i in 1:nT
        S[i,k] ~ TruncatedNormal( (m[i,k] + qc[k]) * q[k], bosd, 1.0e-9, 1.25 )   
      end
    end
    
end


