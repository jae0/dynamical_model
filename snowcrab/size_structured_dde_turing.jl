
using Turing, DifferentialEquations, FillArrays

@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM,
    solver=MethodOfSteps(Tsit5()), dt = 0.01,  ::Type{T} = Float64) where T
  
  # iok= findall(!ismissing, S)
   # biomass process model: 
    K ~ filldist( TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0), nS )  
  
    q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 2.0), nS )    
    qc ~ filldist( TruncatedNormal( 0.0, 0.1, -0.5, 0.5), nS )  

  #  bpsd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.3 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.3 )  ;  # slightly informative .. center of mass between (0,1)
  
    # birth rate from F(y - 8 to 10)  and for males m5 and females
    b ~ filldist( TruncatedNormal(1.0, 0.1, 0.1, 10.0), 2 ) 
    
    # background mortality
    d ~ filldist( TruncatedNormal(0.2, 0.1, 0.1, 0.6), nS )  

    # transition rates
    v ~ filldist( TruncatedNormal(0.8, 0.1, 0.4, 0.99), 4 ) 

    # initial conditions
    u0 ~ filldist( TruncatedNormal( 0.8, 0.1, 0.1, 0.9), nS )

 #    m = Matrix{T}(undef, nM, nS)   # 
   
    pm = ( b, K, d, v, tau, hsa ) 
    # @show pm


    # process model
    msol = solve( 
        remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm ), 
        solver, 
        callback=cb, 
        # maxiters=1e6,
        isoutofdomain=(y,p,t)->any(x->x<0.0, y), 
        saveat=dt
    ) 
    
    # @show msol.retcode
    
    if msol.retcode != :Success 
      Turing.@addlogprob! -Inf
      return nothing
    end
   
    # for i in 1:nM
    #     ii = findall(x->x==prediction_time[i], msol.t)[1]
    #     for k in 1:nS
    #         m[i,k] ~ TruncatedNormal( msol.u[ii][k]/ K[k] , bpsd, 0.01, 0.99 )  
    #     end
  
    # end
  
    for i in 1:nT
        ii = findall(x->x==survey_time[i], msol.t)[1]
        for k in 1:nS
            S[i,k] ~ TruncatedNormal( msol.u[ii][k] / K[k]  * q[k] + qc[k], bosd, 0.0, 1.0)  
        end
    end
     
end
 
