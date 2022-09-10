
using Plots, StatsPlots

function size_structured_dde_turing_plot( ; selection="withfishing withoutfishing S K predictions predictionmeans", si=[1], scale_factor=1.0, mw=nothing )


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
                remake( prob, u0=u0, h=h, tspan=tspan, p=pm; constant_lags=lags ), 
                solver, 
                callback=cb, 
                saveat=dt, 
                isoutofdomain=(y,p,t)->any(x->x<0,y) 
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
            prob2 = DDEProblem( size_structured_dde!, u0, h, tspan, pm; saveat=dt, constant_lags=lags )
            msol = solve( 
                prob2,  
                solver, 
                saveat=dt, 
                isoutofdomain=(y,p,t)->any(x->x<0,y) 
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
        # sample and plot posterior means from model (posterior post-fishery abundance)
        
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
                plot!(prediction_time, w;  alpha=0.01, color=:orange)
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

    plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )

end

