
using Turing, DifferentialEquations

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM, solver=MethodOfSteps(Tsit5())  )     
    
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K ~ filldist( TruncatedNormal( kmu, kmu*0.2, kmu/5.0, kmu*5.0), nS )  
  
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 10.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -0.5, 0.5), nS )  
       
    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.1, 10.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.5, 0.1, 0.25, 0.9), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.5, 1.0), 4 ) 

    # initial conditions
    # starting b prior to first catch event permit higher than 1.0
    m0 ~ filldist( TruncatedNormal(0.9, 0.2, 0.5, 1.25), nS ) 
    u0 = m0 .* K
  
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
     
    msol = solve( prob, 
      solver, 
      callback=cb, 
      saveat=dt,
      isoutofdomain=(u,p,t)->any(x -> x < 0.0, u)
    ) 
     
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    m = tzeros(nM, nS)
    for i in 1:nM
      j = findall(t -> t==prediction_time[i], msol.t)
      if length(j) > 0
        usk =  msol.u[j[1]] ./ K 
        for k in 1:nS
          m[i,k] ~ TruncatedNormal( usk[k], bpsd, 0.0, 1.0)  ; 
        end
      end
    end
  

    # # in case no viable solutions predicted at required time steps
    # if any( x -> x < 0.0, m)
    #   Turing.@addlogprob! -Inf
    #   return nothing
    # end
    
    
    # # check positivity of back transform
    # # column-wise arithmetic
    # yhat = m .+ qc'   
    # if any( x -> x < 0.0, yhat )
    #   Turing.@addlogprob! -Inf
    #   return nothing
    # end
    
    
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for k in 1:nS
      for i in 1:nT
        S[i,k] ~ TruncatedNormal( m[i,k] * q[k] + qc[k] , bosd, 0.0, 1.0 )   
      end
    end
    
    # @show size(m)
    # @show m[28,:]
end


