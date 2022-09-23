
using Turing, DifferentialEquations, FillArrays

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM; 
   jok=ismissing.(S), solver=MethodOfSteps(Tsit5()), dt = 0.01 , ::Type{T} = Float64) where {T}  
  
  # iok= findall(!ismissing, S)
    
    # biomass process model: 
    K ~ filldist( TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0), nS )  
  
    q ~ filldist( TruncatedNormal(  1.0, 0.2,  0.2, 5.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.2, -1.0, 1.0), nS )  
       
    bpsd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.25 )  ;  # slightly informative .. center of mass between (0,1)
  
    # birth rate from F_8 to F_10
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.01, 10.0), 2 ) 
    
    # mortality
    d ~ filldist( TruncatedNormal(0.2, 0.1, 0.01, 1.0), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.9, 0.1, 0.5, 1.0), 4 ) 

    # initial conditions
    # starting b prior to first catch event permit higher than 1.0
    m0 ~ filldist( TruncatedNormal(0.8, 0.2, 0.1, 1.0), nS ) 
    u0 = m0 .* K
  
    p = ( b, K, d, v, tau, hsa )

    # process model
    prob = remake( prob; u0=u0, h=h, tspan=tspan, p=p )

    m = Matrix{T}(undef, nM, nS)
    m = tzeros(nM, nS)  ## predictions

    msol = solve( prob, solver, callback=cb, saveat=dt ) 
     
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end
    
    usk = reduce(hcat, msol.u)'[:,] / K   
    for i in 1:nM
        j = findall(t -> t==prediction_time[i], msol.t)
        if length(j) > 0 
            for k in 1:nS
                m[i,k] ~ TruncatedNormal( usk[j,k][1], bpsd, 0.0, 1.0 )   
            end
        end
    end
    
    # in case no viable solutions predicetd at required time steps
    if any( x -> x < smallnumber, m)
      Turing.@addlogprob! -Inf
      return nothing
    end

    # likelihood at survey times
    # observation model: Y = q X + qc ; X = (Y - qc) / q
    for i in 1:nT
        j = findall(t -> t==survey_time[i], msol.t)
        if length(j) > 0 
        for k in 1:nS
            if jok[i,k] == 0 
                S[i,k] ~ TruncatedNormal( usk[j,k][1] * q[k] + qc[k], bosd, 0.0, 1.0 )   
            end
        end
        end
    end


    # @show size(m)
end



