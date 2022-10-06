
function size_structured_dde!( du, u, h, p, t)
  # here u, du are actual numbers .. not normalized by K due to use of callbacks  

  (b, K, d, v, tau, hsa)  = p

  u1 = h(p, t-1.0)    # no in previous years
  f8 = h(p, t-8.0)[6]  # no mature fem 8  yrs ago
  vh = hsa(t, 1:6)  
  
  # this break down seems to speed it up a bit ... not sure why
  br =  f8 .* b    
  tr =  v .* u1[2:5]
  dr =  d .* u .* u ./ K ./ vh 
  
  du[1] = tr[1]            - dr[1]       # note:       
  du[2] = tr[2]   - tr[1]  - dr[2]     
  du[3] = tr[3]   - tr[2]  - dr[3]    
  du[4] = tr[4]   - tr[3]  - dr[4]    
  du[5] = br[1]   - tr[4]  - dr[5]     
  du[6] = br[2]            - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers
   
end
 



@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM,
  solver=MethodOfSteps(Tsit5()), dt = 0.01,  ::Type{T} = Float64) where T

# iok= findall(!ismissing, S)
 # biomass process model: 
  K ~ filldist( TruncatedNormal( kmu, kmu*0.25, kmu/5.0, kmu*5.0), nS )  

  q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 2.0), nS )    
  qc ~ filldist( TruncatedNormal( 0.0, 0.1, -0.5, 0.5), nS )  

  bosd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.3 )  ;  # slightly informative .. center of mass between (0,1)

  # birth rate from F(y - 8 to 10)  and for males m5 and females
  b ~ filldist( TruncatedNormal(1.0, 0.1, 0.2, 5.0), 2 ) 
  
  # background mortality
  d ~ filldist( TruncatedNormal(0.2, 0.1, 0.05, 0.8), nS )  

  # transition rates
  v ~ filldist( TruncatedNormal(0.8, 0.1, 0.1, 0.99), 4 ) 

  # initial conditions
  u0 ~ filldist( TruncatedNormal( 0.8, 0.1, 0.1, 0.99), nS )

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
  
  for i in 1:nT
      ii = findall(x->x==survey_time[i], msol.t)[1]
      for k in 1:nS
          S[i,k] ~ TruncatedNormal( msol.u[ii][k] / K[k]  * q[k] + qc[k], bosd, 0.0, 1.0)  
      end
  end
   
end




function size_structured_predictions_annual( res; prediction_time=prediction_time, n=1e10 )
 
    nchains = size(res)[3]
    nsims = size(res)[1]
    
    nI = min( nchains*nsims , n )
  
    m = zeros(nM, nS, nI, 2)
   
    z = 0
    
    for j in 1:nsims  # nsims 
    for l in 1:nchains #nchains
  
      z += 1
  
      b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
      K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
      d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
      v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
  
      q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
      qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]
  
      u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]
  
      pm = ( b, K, d, v, tau, hsa ) 
      
      prb = remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm ) 
      
      msol = solve( prb, solver, callback=cb, saveat=dt ) 
      msol2 = solve( prb, solver, saveat=dt ) # no call backs
  
      for i in 1:nM
          ii = findall(x->x==prediction_time[i], msol.t)[1]
          jj = findall(x->x==prediction_time[i], msol2.t)[1]
          m[i,:,z,1] = msol.u[ii]    # with fishing   
          m[i,:,z,2] = msol2.u[jj]   # no fishing
      end
  
      z >= nI && return m         
  
    end
    end
  
    return m 
    
end    


# -----------


function size_structured_predictions( res; n=10, k=1 )
 
  nchains = size(res)[3]
  nsims = size(res)[1]

  z = 0

  gr() 
  theme(:default)
  pl =plot()
  
  for j in 1:nsims  # nsims 
  for l in 1:nchains #nchains
      z += 1

      b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
      K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
      d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
      v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]

      q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
      qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]

      u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

      pm = ( b, K, d, v, tau, hsa ) 
      
      prb = remake( prob; u0=u0 .* K , h=h, tspan=tspan, p=pm ) 

      msol2 = solve( prb, solver, saveat=dt ) # no call backs
      
      sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t) ./ 1000.0 ./ 1000.0 :  scale_factor

      yval2 = vec( reduce(hcat, msol2.u)'[:,k] .* sf2)
      pl = plot!( pl, msol2.t, yval2, alpha=0.01, lw=1, color=:lime ) 

      if k==1
          # do fishing 
          msol = solve( prb, solver, callback=cb, saveat=dt ) 

          sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t) ./ 1000.0 ./ 1000.0 :  scale_factor
  
          yval = vec( reduce(hcat, msol.u)'[:,k] .* sf )
          pl = plot!( pl, msol.t, yval, alpha=0.01, lw=1, color=:orange ) 
      end
      
      if z >= n 
          pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
          # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
          pl =  plot!(pl; legend=false )

          return pl
      end

  end
  end

  pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
  pl =  plot!(pl; legend=false )

  return pl
end    
   

function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end
 

function removals_aggregate( removed, fish_time )
  f = DataFrame( yr=floor.(fish_time), rem = removed );
  out = combine(groupby(f,:yr),[:rem ] .=> sum )
  return(out)
end


function fishing_mortality( removed, fish_time, biomass  )
  removed_annual = removals_aggregate( removed, fish_time )
  Fkt = removed_annual[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  FR =  Fkt ./ ( Fkt .+  biomass )  # relative F
  FM = -log.(  1.0 .- ( FR ) )  # instantaneous F
  return ( Fkt, FR, FM )
end



function size_structured_dde_turing_plot( ; selection="withfishing withoutfishing S K predictions predictionmeans", si=[1], scale_factor=1.0, mw=nothing, vn="" )
  
   
  if occursin( r"withfishing", selection ) | occursin( r"withoutfishing", selection )  
      # mean field dynamics:
      u0 = [ 
          mean( res[:,"K[1]",:] ),
          mean( res[:,"K[2]",:] ),
          mean( res[:,"K[3]",:] ),
          mean( res[:,"K[4]",:] ),
          mean( res[:,"K[5]",:] ),
          mean( res[:,"K[6]",:] )
          ] .*  [ 
          mean( res[:,"m[1,1]",:] ),
          mean( res[:,"m[1,2]",:] ),
          mean( res[:,"m[1,3]",:] ),
          mean( res[:,"m[1,4]",:] ),
          mean( res[:,"m[1,5]",:] ),
          mean( res[:,"m[1,6]",:] )
      ]

      b = [ mean( res[:,"b[1]",:] ), mean( res[:,"b[2]",:] ) ]
      K = [ mean( res[:,"K[1]",:] ), mean( res[:,"K[2]",:] ), 
            mean( res[:,"K[2]",:] ), mean( res[:,"K[2]",:] ),
            mean( res[:,"K[2]",:] ), mean( res[:,"K[2]",:] ) ]   ; 
      d = [ mean( res[:,"d[1]",:] ), mean( res[:,"d[2]",:] ), 
          mean( res[:,"d[3]",:] ), mean( res[:,"d[4]",:] ),
          mean( res[:,"d[5]",:] ), mean( res[:,"d[6]",:] ) ]   
      v = [ mean( res[:,"v[1]",:] ), mean( res[:,"v[2]",:] ), 
          mean( res[:,"v[3]",:] ), mean( res[:,"v[4]",:] ) ]

      pm = ( b, K, d, v, tau, hsa ) 
      
      if occursin( r"withfishing", selection )

          msol = solve( 
              remake( prob, u0=u0, h=h, tspan=tspan, p=pm; constant_lags=[tau] ), 
              solver, 
              callback=cb, 
              saveat=dt #, 
              # isoutofdomain=(y,p,t)->any(x->x<0,y) 
          )
          yval = reduce(hcat, msol.u)'[:,si]

          if nameof(typeof(mw)) == :ScaledInterpolation
              yval = yval .* mw(msol.t) ./ 1000.0  ./ 1000.0 
          else
              yval = yval .* scale_factor
          end
          yval = vec(yval)
          
          plot!( msol.t, yval, alpha=0.85, lw=4, color=:steelblue ) 
          plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
          plot!(; ylim=(0, maximum(yval)*1.1 ) )
          plot!(; legend=false )
      end

  
      if occursin( r"withoutfishing", selection )  
          prob2 = DDEProblem( size_structured_dde!, u0, h, tspan, pm; saveat=dt, constant_lags=[tau] )
          msol = solve( 
              prob2,  
              solver, 
              saveat=dt#, 
              # isoutofdomain=(y,p,t)->any(x->x<0,y) 
          )
          yval = reduce(hcat, msol.u)'[:,si]

          if nameof(typeof(mw)) == :ScaledInterpolation
              yval = yval .* mw(msol.t) ./ 1000.0  ./ 1000.0 
          else
              yval = yval .* scale_factor
          end
          yval = vec(yval)

          plot!( msol.t, yval, alpha=0.85, lw=4, color=:teal ) 
          plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
          plot!(; ylim=(0, maximum(yval)*1.1 ) )
          plot!(; legend=false )
      end
  end

  if occursin( r"K", selection )  
      # sample and plot posterior K
      si1 = si[1]

      if nameof(typeof(mw)) == :ScaledInterpolation
          sf = mean(mw(yrs)) ./ 1000.0  ./ 1000.0 
          for i in 1:length(res)  
              w = res[i,Symbol("K[$si1]"),1]  .* sf
              hline!([w];  alpha=0.05, color=:limegreen )
          end
          o = vec( reduce(hcat, res[:,Symbol("K[$si1]"),:]) ) .* sf
      else
          for i in 1:length(res)  
              w = res[i,Symbol("K[$si1]"),1]  .* scale_factor
              hline!([w];  alpha=0.05, color=:limegreen )
          end
          o = vec( reduce(hcat, res[:,Symbol("K[$si1]"),:]) ) .* scale_factor

      end
      
      hline!([mean(o)];  alpha=0.6, color=:chartreuse4, lw=5 )
      hline!([quantile(o, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
      hline!([quantile(o, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

      plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
      plot!(; ylim=(0, maximum(o)*1.1 ) )

      plot!(; legend=false )
  end


  if occursin( r"predictions", selection )  
      # sample and plot posterior means from model (posterior abundance   at sept )
      
      for k in si
          for l in 1:size(res)[3]
          # l = 1
          for i in 1:length(res)  
              w = zeros(nM)
              for j in 1:nM
                  w[j] = res[i, Symbol("K[$k]"),l] * res[i, Symbol("m[$j,$k]"),l] 
              end
              if nameof(typeof(mw)) == :ScaledInterpolation
                  w = w .* mw(prediction_time) ./ 1000.0  ./ 1000.0 
              else
                  w = w .* scale_factor
              end
              plot!(prediction_time, w;  alpha=0.1, color=:orange)
          end
          end
      end
      plot!(; legend=false )
  end

  
  if occursin( r"predictionmeans", selection )  
      # mean post-fishery abundance
      
      for k in si
          u = zeros(nM)
          v = zeros(nM)
          for  j in 1:nM
              u[j] = mean( res[:,Symbol("m[$j,$k]"),:] .* res[:,Symbol("K[$k]"),:] )   
              v[j] = std(  res[:,Symbol("m[$j,$k]"),:] .* res[:,Symbol("K[$k]"),:] )  
          end
          if nameof(typeof(mw)) == :ScaledInterpolation
              u = u .* mw(prediction_time) ./ 1000.0  ./ 1000.0 
              v = v .* mw(prediction_time) ./ 1000.0  ./ 1000.0 
          else
              u = u .* scale_factor
              v = v .* scale_factor
          end
          plot!(prediction_time, u, lw=2, color=:orangered )
          scatter!(prediction_time, u, markersize=4, color=:goldenrod1 )
          
      end
  end


  if occursin( r"S", selection )  
      # back transform S to normal scale 
      si1 = si[1]
      yhat = ( S[:,si1] ./ mean(res[:,Symbol("q[$si1]"),:]) .- mean(res[:,Symbol("qc[$si1]"),:] ) ) .* mean(res[:,Symbol("K[$si1]"),:]  ) 
      if nameof(typeof(mw)) == :ScaledInterpolation
          yhat = yhat .* mw(yrs) ./ 1000.0  ./ 1000.0 
      else
          yhat = yhat .* scale_factor
      end
      plot!(survey_time, yhat, color=:purple2, lw=2 )
      scatter!(survey_time, yhat, markersize=4, color=:purple4)
      plot!(; legend=false )
  end


end


function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end
  