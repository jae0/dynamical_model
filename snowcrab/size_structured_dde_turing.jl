
using Turing, DifferentialEquations

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM; 
   jok=ismissing.(S), solver=MethodOfSteps(Tsit5()), saveat = 0.01  )     
  
  # iok= findall(!ismissing, S)
    
    # biomass process model: 
    K ~ filldist( TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0), nS )  
  
    q ~ filldist( TruncatedNormal(  1.0, 0.25,  0.2, 5.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.25, -1.0, 1.0), nS )  
       
    bpsd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.25, 0.1, 10.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 0.5), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.5, 1.0), 4 ) 

    # initial conditions
    # starting b prior to first catch event permit higher than 1.0
    m0 ~ filldist( TruncatedNormal(0.9, 0.2, 0.5, 1.0), nS ) 
    u0 = m0 .* K
  
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
  
    m = tzeros(nM, nS)
    ms = tzeros(nT, nS) ## survey estimates
    
    msol = solve( prob, 
      solver, 
      callback=cb, 
      # isoutofdomain=(u,p,t)->any(x -> x < smallnumber, u),  # do not use as it seems to cause issues with stability of turing
      saveat=saveat
      # ,
      # tstops=survey_time  # make sure solutions exist at prediction_time
    ) 
     
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
    
    for i in 1:nM
      j = findall(t -> t==prediction_time[i], msol.t)
      m[i,:] ~ MvNormal( msol.u[j[1]] ./ K, bpsd )
    end

    # likelihood at survey times
    for i in 1:nT
      j = findall(t -> t==survey_time[i], msol.t)
      ms[i,:] ~ MvNormal( msol.u[j[1]] ./ K, bpsd )
    end

    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for k in 1:nS
      for i in 1:nT
        if jok[i,k] == 0 
            S[i,k] ~ TruncatedNormal( ms[i,k] * q[k] + qc[k], bosd, smallnumber, 1.0 )   
        end
      end
    end
    
    # @show size(m)
end


