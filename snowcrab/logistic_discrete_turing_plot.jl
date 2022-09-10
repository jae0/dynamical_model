
using Plots, StatsPlots

function logistic_discrete_turing_plot( ; selection="S K predictions predictionmeans"  )

 
    if occursin( r"K", selection )  
        # sample and plot posterior K

        for i in 1:length(res)  
            w = res[i,Symbol("K"),1]   
            hline!( [w];  alpha=0.05, color=:limegreen )
        end
        o = vec( reduce(hcat, res[:,Symbol("K"),:]) )  
        
        hline!([mean(o)];  alpha=0.6, color=:chartreuse4, lw=5 )
        hline!([quantile(o, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
        hline!([quantile(o, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )

        plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+1.5  ) )
        plot!(; ylim=(0, maximum(o)*1.1 ) )

        plot!(; legend=false )
    end


    if occursin( r"predictions", selection )  
        # sample and plot posterior means from model (posterior post-fishery abundance)

        for l in 1:size(res)[3]
            for i in 1:length(res)  
                w = zeros(nM)
                for j in 1:nM
                    w[j] = res[i, Symbol("K"),l] * res[i, Symbol("m[$j]"),l] 
                end

                plot!(prediction_time, w;  alpha=0.01, color=:orange)
            end
        end

        plot!(; legend=false )
    end

    
    if occursin( r"predictionmeans", selection )  
        # mean post-fishery abundance

        u = zeros(nM)
        v = zeros(nM)
        for  j in 1:nM
            u[j] = mean( res[:,Symbol("m[$j]"),:] .* res[:,Symbol("K"),:] )   
            v[j] = std(  res[:,Symbol("m[$j]"),:] .* res[:,Symbol("K"),:] )  
        end
        plot!(prediction_time, u, lw=2, color=:orangered )
        scatter!(prediction_time, u, markersize=4, color=:goldenrod1 )
    end


    if occursin( r"S", selection )  
        # back transform S to normal scale 
        if model_variation == "Model_2" 
            # observation model: S[Y] ~ m[X] --> Y = q X + qc ; 
            # inverse operation: X = (Y - qc) / q
            yhat = ( S[:] .- mean(res[:,Symbol("qc"),:] ) ) ./ mean(res[:,Symbol("q"),:]) .* mean(res[:,Symbol("K"),:]  ) 
        else
            # observation model: Y = q X  ; 
            # inverse operation: X = (Y ) / q
            yhat = ( S[:] ./ mean(res[:,Symbol("q"),:]) ) .* mean(res[:,Symbol("K"),:]  ) 
        end

        plot!(survey_time, yhat, color=:purple2, lw=2 )
        scatter!(survey_time, yhat, markersize=4, color=:purple4)
        plot!(; legend=false )

    end

    plot!(; xlim=(minimum(yrs)-0.5, maximum(yrs)+nP  ) )

end


  