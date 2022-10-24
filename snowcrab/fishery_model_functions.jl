
# generic functions shared by fishery_models

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


function plots_diagnostic( res, vn="K" ) 
  gr()
  pl = plot()
  pl = density!(pl, res[ Symbol(vn) ])
  return pl
  # vn = "b[1]"; density!(res[ Symbol(vn) ])
  # vn = "b[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[1]"; density!(res[ Symbol(vn) ])
  # vn = "d[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[3]"; density!(res[ Symbol(vn) ])
  # vn = "d[4]"; density!(res[ Symbol(vn) ])
  # vn = "d[5]"; density!(res[ Symbol(vn) ])
  # vn = "d[6]"; density!(res[ Symbol(vn) ])
  # vn = "K[1]"; density!(res[ Symbol(vn) ])
  # vn = "K[2]"; density!(res[ Symbol(vn) ])
  # vn = "K[3]"; density!(res[ Symbol(vn) ])
  # vn = "K[4]"; density!(res[ Symbol(vn) ])
  # vn = "K[5]"; density!(res[ Symbol(vn) ])
  # vn = "K[6]"; density!(res[ Symbol(vn) ])
end



# ----------


function fishery_model_inference( fmod; n_adapts=1000, n_samples=1000, n_chains=4, max_depth=7, init_ϵ=0.05, 
  turing_sampler = Turing.NUTS(n_adapts, 0.65; max_depth=max_depth, init_ϵ=init_ϵ), debug=false, seed=1  )
  
  if !debug 
    Logging.disable_logging(Logging.Warn) # or e.g. Logging.Info
  else
    n_chains=1  
  end
  
  Random.seed!(seed)
  
  # 1000 -> ? hrs (Tsit5);  500 -> 6 hrs;; 29hrs 100/100 cfasouth
  #   # n_chains = Threads.nthreads()
  # turing_sampler = Turing.NUTS(n_adapts, 0.65; max_depth=max_depth, init_ϵ=init_ϵ)  ;# stepsize based upon previous experience
  
  res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains )
  # if on windows and threads are not working, use single processor mode:
  # res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

  showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

  return res
end



# -----------


function fishery_model_mortality( removed, fb )    
  if @isdefined fish_year
    removed_annual_kt = removals_aggregate( removed, fish_year )
    Fkt = removed_annual_kt[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  else 
    Fkt = removed
  end
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  pl = plot()
  pl = plot!(pl, survey_time, FM ;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, survey_time, mean(FM, dims=2) ;  alpha=0.8, color=:slateblue, lw=4)
  pl = plot!(pl, xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
  pl = plot!(pl, ylim=(0, maximum(FM)*1.1 ) )
  pl = plot!(pl ; legend=false )
  return ( Fkt, FR, FM, pl )
end



# -----------


function fishery_model_harvest_control_rule(res, yrs; FM=FM, fb=fb, n_sample=500 )

  pl = plot()

  if  occursin( r"size_structured", model_variation ) 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor

    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass
  else
    K = vec( Array(res[:, Symbol("K"), :]) ) 
  end
  
  o = rand(K, n_sample)
  pl = vline!(pl, o;  alpha=0.05, color=:limegreen )
  pl = vline!(pl, o./2;  alpha=0.05, color=:darkkhaki )
  pl = vline!(pl, o./4;  alpha=0.05, color=:darkred )

  pl = vline!(pl, [mean(o)];  alpha=0.6, color=:chartreuse4, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

  pl = vline!(pl, [mean(o)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )

  pl = vline!(pl, [mean(o)/4.0];  alpha=0.6, color=:darkred, lw=5 )
  pl = vline!(pl, [quantile(o, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  pl = vline!(pl, [quantile(o, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )

  nt = length(survey_time)
  colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]

  # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)

  fb_mean = mean(fb, dims=2)
  fm_mean = mean(FM, dims=2)

  # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
  pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)

  pl = scatter!(pl,  fb_mean, fm_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
    series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
  pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
  ub = max( quantile(o, 0.95), maximum( fb_mean ) ) * 1.05
  pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
  # TODO # add predictions ???
  return(K, fb_mean, fm_mean, pl)
end

