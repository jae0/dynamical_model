

function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks

    (b, K, d, v, a, tau, hsa)  = p

    u1 = h(p, t-1.0)[2:5]    # no in previous years
    f8 = h(p, t-8.0)[6]  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6)

    # this break down seems to speed it up a bit ... not sure why
    br =  f8 .* b
    tr =  v .* u1
    # dr =  d .* ( max.(u,0.0) ./ vh) .^ (2.0)

    uv = u ./ vh
    dr =  d .* uv .* uv

    # a = K[2:6] ./ K[1:5]

    du[1] = tr[1] * a[1]            - dr[1]       # note:
    du[2] = tr[2] * a[2]   - tr[1]  - dr[2]
    du[3] = tr[3] * a[3]   - tr[2]  - dr[3]
    du[4] = tr[4] * a[4]   - tr[3]  - dr[4]
    du[5] = br[1] * a[5]   - tr[4]  - dr[5]
    du[6] = br[2]                   - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers

end


function dde_parameters()
   # these are dummy initial values .. just to get things started
  b=[1.7, 0.5]
  K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
  a=[1.0, 1.0, 1.0, 1.0, 1.0];

  d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];

  v=[0.65, 0.68, 0.61, 0.79];

   params = ( b, K, d, v, a,  tau, hsa)
   return params
end


@model function size_structured_dde_turing( S, kmu, tspan, prob, nT, nS, nM,
  solver=MethodOfSteps(Tsit5()), dt = 0.01,  ::Type{T} = Float64) where T

  # biomass process model:
  K ~ filldist( TruncatedNormal( kmu, kmu*0.2, kmu/10.0, kmu*5.0), nS )  # kmu is max of a multiyear group , serves as upper bound for all

  a ~ MvNormal( K[2:6] ./ K[1:5], 0.5 )

  q ~ filldist( TruncatedNormal(  1.0, 0.1,  0.1, 2.0), nS )
  qc ~ filldist( TruncatedNormal( 0.0, 0.1, -0.5, 0.5), nS )

  model_sd ~  TruncatedNormal( 0.1, 0.1, 0.01, 0.35)

  # birth rate from F(y - 8 to 10)  and for males m5 and females
  b ~ filldist( TruncatedNormal(1.0, 0.2, 0.01, 100.0), 2 )

  # background mortality
  d ~ filldist( TruncatedNormal(0.4, 0.1, 0.01, 0.99), nS )

  # transition rates
  v ~ filldist( TruncatedNormal(0.8, 0.1, 0.01, 0.99), 4 )

  # initial conditions
  u0 ~ filldist( TruncatedNormal(0.5, 0.1, 0.01, 0.99), nS )

  pm = ( b, K, d, v, a, tau, hsa )
  # @show pm

  # m = TArray{T}( length(Si), nS )

  # process model
  msol = solve(
      remake( prob; u0=u0, h=h, tspan=tspan, p=pm ),
      solver,
      callback=cb,
      # maxiters=1e6,
      isoutofdomain=(y,p,t)->any(x -> (x<0.0 || x>1.0), y) ,
      #isoutofdomain=(y,p,t)->any(x -> (x<0.0), y),
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    # fill this as without causes Turing confusion
    # for i in 1:nSI
    #     for k in 1:nS
    #         m[i,k] ~ TruncatedNormal( smallnumber, smallnumber , 0.0, 1.0 )  # dummy values as m has not been sampled yet .. without this it will crash
    #     end
    # end
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ TruncatedNormal( msol.u[ii][k] * q[k] + qc[k], model_sd, 0.0, 1.0)  # observation and process error combined
      end
  end

end




function size_structured_predictions_annual( res; prediction_time=prediction_time, ns=1e10 )

    nchains = size(res)[3]
    nsims = size(res)[1]

    nI = Int( min( nchains*nsims , ns ) )

    md = zeros(nM, nS, nI, 2)  # number normalized
    mn = zeros(nM, nS, nI, 2)  # numbers
    mb = mn[:,1,:,:]  # biomass of first class

    z = 0

    for j in 1:nsims  # nsims
    for l in 1:nchains #nchains

      z += 1

      b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
      K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
      a = [ res[j, Symbol("a[$k]"), l] for k in 1:(nS-1)]

      d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]

      v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]

      q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
      qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]

      u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

      pm = ( b, K, d, v, a, tau, hsa )

      prb = remake( prob; u0=u0 , h=h, tspan=tspan, p=pm )

      msol = solve( prb, solver, callback=cb, saveat=dt )
      msol2 = solve( prb, solver, saveat=dt ) # no call backs

      for i in 1:nM
          ii = findall(x->x==prediction_time[i], msol.t)[1]
          jj = findall(x->x==prediction_time[i], msol2.t)[1]
          sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t[ii])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
          sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t[jj]) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
          md[i,:,z,1] = msol.u[ii]   # with fishing
          md[i,:,z,2] = msol2.u[jj]  # no fishing
          mn[i,:,z,1] = msol.u[ii]   .* K  # with fishing
          mn[i,:,z,2] = msol2.u[jj]   .* K # no fishing
          mb[i,z,1] = mn[i,1,z,1]  .* sf
          mb[i,z,2] = mn[i,1,z,2]  .* sf2

      end

      z >= nI && return (md, mn, mb)

    end
    end

    return (md, mn, mb)

end


# -----------


function size_structured_predictions( res; ns=10, plot_k=1, alpha=0.01 )

  nchains = size(res)[3]
  nsims = size(res)[1]

  nI = Int( min( nchains*nsims , ns ) )

  z = 0

  out = Vector{Vector{Float64}}()
  out2 = Vector{Vector{Float64}}()

  gr()
  theme(:default)
  pl =plot()

  for j in 1:nsims  # nsims
  for l in 1:nchains #nchains
      z += 1

      b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
      K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
      a = [ res[j, Symbol("a[$k]"), l] for k in 1:(nS-1)]

      d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]

      v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]

      q =  [ res[j, Symbol("q[$k]"), l] for k in 1:nS]
      qc = [ res[j, Symbol("qc[$k]"), l] for k in 1:nS]

      u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

      pm = ( b, K, d, v, a, tau, hsa )

      prb = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

      if plot_k==1
          # do fishing and nonfishing

          msol = solve( prb, solver, callback=cb, saveat=dt )
          msol2 = solve( prb, solver, saveat=dt ) # no call backs

          sf = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol.t) ./ 1000.0 ./ 1000.0 :  scale_factor
          sf2 = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol2.t) ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt

          yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k]) .* K[plot_k] .* sf2
          yval = vec( reduce(hcat, msol.u)'[:,plot_k] ) .* K[plot_k] .* sf

          pl = plot!( pl, msol.t, yval, alpha=alpha, lw=1, color=:orange )
          pl = plot!( pl, msol2.t, yval2, alpha=alpha*2.0, lw=1, color=:lime )

          push!(out, yval)
          push!(out2, yval2)

      else
          msol2 = solve( prb, solver, saveat=dt ) # no call backs
          yval2 = vec( reduce(hcat, msol2.u)'[:,plot_k]) .* K[plot_k]
          pl = plot!( pl, msol2.t, yval2, alpha=alpha, lw=1, color=:lime )

          push!(out2, yval2)
      end

      if z >= nI
          pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
          # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
          pl =  plot!(pl; legend=false )

          return (out, out2, pl)
      end

  end
  end

  pl =  plot!(pl; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  # pl =  plot!(pl; ylim=(0, maximum(m[:,:,2,z])*1.1 ) )
  pl =  plot!(pl; legend=false )

  return (out, out2, pl)
end


function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end


function removals_aggregate( removed, fish_year )
  landings_aggregated = DataFrame( yr=floor.(fish_year), rem = removed );
  out = combine(groupby(landings_aggregated,:yr),[:rem ] .=> sum )
  stimes = DataFrame( yr=floor.(survey_time) )
  out = leftjoin(stimes, out, on=:yr)
  sort!(out, :yr)
  oo = findall(x->ismissing(x), out[:,:rem_sum])
  if length(oo) > 0
    out[ oo[1], :rem_sum ] = 0.0
  end
  return(out)
end



function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end

