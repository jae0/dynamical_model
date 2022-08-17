

function plot_dde( ; selection="withfishing withoutfishing S K predictions predictionmeans", s=[1] )

    if occursin( r"withfishing", selection )  
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
            mean( res[:,"K[2]",:] ), mean( res[:,"K[2]",:] ) ]  ; 
        d = [ mean( res[:,"d[1]",:] ), mean( res[:,"d[2]",:] ), 
            mean( res[:,"d[3]",:] ), mean( res[:,"d[4]",:] ),
            mean( res[:,"d[5]",:] ), mean( res[:,"d[6]",:] ) ]   
        v = [ mean( res[:,"v[1]",:] ), mean( res[:,"v[2]",:] ), 
            mean( res[:,"v[3]",:] ), mean( res[:,"v[4]",:] ) ]

        pm = ( b, K, d, v, tau, hsa ) 

        msol = solve( 
            remake( prob, u0=u0, h=h, tspan=tspan, p=pm ), 
            solver, 
            callback=cb, 
            saveat=dt, 
            isoutofdomain=(y,p,t)->any(x->x<0,y) 
        )
        plot!( msol.t, reduce(hcat, msol.u)'[:,s], alpha=0.5, lw=5 ) 
        plot!(; legend=false )
    end

    if occursin( r"withoutfishing", selection )  
        prob2 = DDEProblem( size_structured!, u0, h, tspan, pm, saveat=dt )
        msol = solve( 
            prob2,  
            solver, 
            saveat=dt, 
            isoutofdomain=(y,p,t)->any(x->x<0,y) 
        )
        plot!( msol.t, reduce(hcat, msol.u)'[:,s], alpha=0.75, lw=5 ) 
        plot!(; legend=false )
    end

 
    if occursin( r"S", selection )  
        # back transform S to normal scale 
        s1 = 1
        yhat = ( S[:,s1] .* mean(res[:,Symbol("q[$s1]"),:]) .- mean(res[:,Symbol("qc[$s1]"),:] ) ) .* mean(res[:,Symbol("K[$s1]"),:] ) 
        scatter!(survey_time, yhat, markersize=4)
        plot!(survey_time, yhat )
        plot!(; legend=false )
    end

    if occursin( r"K", selection )  
        # sample and plot posterior K
        s1 = 1  # state variable index
        for i in 1:length(res)  
            w = res[i,Symbol("K[$s1]"),1]
            hline!([w];  alpha=0.1 )
        end
        plot!(; legend=false )
    end


    if occursin( r"predictions", selection )  
        # sample and plot posterior means from model (posterior post-fishery abundance)
        for k in s
            w = zeros(nT)
            for i in 1:length(res)  
                for j in 1:nT
                    w[j] = res[i, Symbol("K[$k]"),1] * res[i, Symbol("m[$j,$k]"),1]
                end
                plot!(prediction_time, w;  alpha=0.1)
            end
        end
        plot!(; legend=false )
    end

    
    if occursin( r"predictionmeans", selection )  
        # mean post-fishery abundance
        for k in s
            u = zeros(nT)
            v = zeros(nT)
            for  j in 1:nT
                u[j] = mean( res[:,Symbol("m[$j,$k]"),:] .* res[:,Symbol("K[$k]"),:] ) 
                v[j] = std(  res[:,Symbol("m[$j,$k]"),:] .* res[:,Symbol("K[$k]"),:] ) 
            end
            scatter!(prediction_time, u, markersize=2 )
        end
        plot!(; xlim=(minimum(yrs), maximum(yrs)  ) )
    end
end
