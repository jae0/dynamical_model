
using Turing, DifferentialEquations

@model function fishery_model_turing_dde( S, kmu, tspan, prob, nT, nS, solver=MethodOfSteps(Tsit5()), force_positive=true )
    
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K ~ filldist( TruncatedNormal( kmu, kmu*0.1, kmu/5.0, kmu*5.0), nS )  
    
    # consider more diffuse Cauchy prior for k .. slow mixing
    # K ~ filldist( truncated( Cauchy( kmu, kmu*0.1), kmu/4.0, kmu*4.0), nS )  
    
    bpsd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.05, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.2, 5.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -0.5, 0.5), nS )  
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.2, 5.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.5, 0.2, 0.1, 0.9), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.8, 0.2, 0.4, 1.0), 4 ) 

    # initial conditions
    m = TArray{Float64}(nT, nS)
    for k in 1:nS 
      m[1,k] ~  TruncatedNormal( 0.8, 0.2, 0.4, 1.25 )  ; # starting b prior to first catch event permit higher than 1.0
    end 

    u0 = [ m[1,1], m[1,2], m[1,3], m[1,4], m[1,5], m[1,6]  ] .* K  # don't know why but takiing an array slice of m causes an error
     # u0 = [vec(m[1,:])] .* K  # don't know why but taking an array slice of m causes an error
    
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
     
    msol = solve( prob, 
      solver, 
      callback=cb, 
      saveat=dt, 
      abstol=1e-9,
      reltol=1e-9,
      isoutofdomain=(u,p,t)->any(x -> x < eps, u)  
    ) 
     
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
 
    # process error model
    for i in 2:nT
      j = findall(t -> t==prediction_time[i], msol.t)
      if length(j) > 0
        usk = max.( msol.u[j[1]], 1.0) ./ K 
        for k in 1:nS
          try
            m[i,k] ~ TruncatedNormal( usk[k], bpsd, 1.0e-9, 1.25)  ; 
          catch
            Turing.@addlogprob! -Inf
            return nothing
          end    
        end
      end
    end

    # observation error model transform to scale of S
  
     for k in 1:nS
      for i in 1:nT
        o = (m[i,k] + qc[k]) * q[k]
        try
          S[i,k] ~ TruncatedNormal( o, bosd, 1.0e-9, 1.25 )   
        catch
          Turing.@addlogprob! -Inf
          return nothing
        end    
      end
    end
    
end

