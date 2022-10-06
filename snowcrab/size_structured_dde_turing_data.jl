# perpare dat for dde run of fishery model
using DifferentialEquations, Interpolations, RData


fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
o = load( fndat, convert=true)
Y = o["Y"]
Kmu = o["Kmu"]
removals = o["L"]
MW = o["M0_W"]
 

    if false
        # alternatively, if running manually:
        # can run R-code that creates local RData file with required data
        # run in R externally or from within julia or .. 
        # from within julia
    
        # using RCall
        # # type $ in Julia's command prompt starts an R session.  
        # # .. run below
        # # type <backspace> to escape back to julia
        
        # source( file.path( code_root, "bio_startup.R" )  )
        # require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        # fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics",  for_julia=TRUE, time_resolution=1/12)
        
        # # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
        # @rget Y  
        # @rget Kmu 
        # @rget removals 
        # @rget ty

        # example line plots
        plot(  Y[:,:yrs], Y[:,:cfasouth_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfasouth_M1] )
        plot!( Y[:,:yrs] .+2, Y[:,:cfasouth_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfasouth_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfasouth_M4] )

        plot(  Y[:,:yrs], Y[:,:cfanorth_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfanorth_M1] )
        plot!( Y[:,:yrs] .+2 , Y[:,:cfanorth_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfanorth_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfanorth_M4] )


        plot(  Y[:,:yrs], Y[:,:cfa4x_M0] )
        plot!( Y[:,:yrs] .+1 , Y[:,:cfa4x_M1] )
        plot!( Y[:,:yrs] .+2 , Y[:,:cfa4x_M2] )
        plot!( Y[:,:yrs] .+3, Y[:,:cfa4x_M3] )
        plot!( Y[:,:yrs] .+4, Y[:,:cfa4x_M4] )
    end
 


# "survey index"
statevars = [
  Symbol("$aulab","_M0"),
  Symbol("$aulab","_M1"),
  Symbol("$aulab","_M2"),
  Symbol("$aulab","_M3"),
  Symbol("$aulab","_M4"),
  Symbol("$aulab","_f_mat")
]

S = Matrix(Y[:, statevars ])
   
# interpolating function for mean weight
mwspline = extrapolate( interpolate( MW[:,Symbol("mw_", "$aulab") ], (BSpline(Linear()) ) ),  Interpolations.Flat() )
mw = Interpolations.scale(mwspline, yrs )

scale_factor = mw(yrs) / (1000 *1000 ) # convert numbers to kt biomass , also used in plots 


# convert to (biomass kt to number) 

# id index
ki = aulab=="cfanorth" ? 1 :
     aulab=="cfasouth" ? 2 :
     aulab=="cfa4x"    ? 3 :
     0  # default

kmu  =  Kmu[ki] / mean(scale_factor)

smallnumber = 1.0 / kmu / 100.0  # floating point value of sufficient to assume 0 valued
    
no_digits = 3  # time floating point rounding 

# spin up time of ~ 1 cycle prior to start of dymamics and project nP years into the future
tspan = (minimum(yrs) - 8.1, maximum(yrs) + nP + 1.1 )  

survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)    # time of observations for survey
  
# this only adds habitat space  ... predation is also a useful one .. 
# speed is the issue 
 
prediction_time = 
  floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  #yrs
  round(round( 9.0/12.0 /dt; digits=0 ) *dt; digits=no_digits)   # sept

#  sa to fraction
external_forcing =  reshape( [
   Y[:,Symbol("H", "$aulab","_M0")]  / maximum( Y[:,Symbol("H", "$aulab","_M0")] )
   Y[:,Symbol("H", "$aulab","_M1")]  / maximum( Y[:,Symbol("H", "$aulab","_M1")] )
   Y[:,Symbol("H", "$aulab","_M2")]  / maximum( Y[:,Symbol("H", "$aulab","_M2")] )
   Y[:,Symbol("H", "$aulab","_M3")]  / maximum( Y[:,Symbol("H", "$aulab","_M3")] )
   Y[:,Symbol("H", "$aulab","_M4")]  / maximum( Y[:,Symbol("H", "$aulab","_M4")] )
   Y[:,Symbol("H", "$aulab","_f_mat")]  / maximum( Y[:,Symbol("H", "$aulab","_f_mat")] ) 
  ], nT, nS )


efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
hsa = Interpolations.scale(efc, yrs .+ 0.75, 1:nS )

fish_time =  round.( round.( removals[:,:ts] ./ dt; digits=0 ) .* dt; digits=no_digits)    # time of observations for survey
removed = removals[:,Symbol("$aulab")]

# keep nonzero elements
i = findall( x-> x>0, removed)
removed = removed[ i ]
fish_time = fish_time[i]


function affect_fishing!(integrator)
  i = findall(t -> t == integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
  # / integrator.p[2][1]  # p[2] ==K divide by K  .. keep unscaled to estimate magnitude of other components
end

# cb = CallbackSet( 
#   PresetTimeCallback( fish_time, affect_fishing! ), 
#   PositiveDomain()
# );

cb = PresetTimeCallback( fish_time, affect_fishing! ) 

# history function 0.5 default
# h(p,t) = ones( nS ) .* 0.5  #values of u before t0
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* kmu
  
tau = 1.0  # delay resolution
 
# these are dummy initial values .. just to get things started

u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .*kmu; 
b=[1.0, 0.8]
K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu; 
d=[0.2, 0.3, 0.4, 0.5, 0.5, 0.5];

v=[0.8, 1.0, 1.0, 1.0];   
 

# dummy values needed to start the turing initialization (bootstrap it)
p = ( b, K, d, v, tau, hsa)    

if false
  ## test, can ignore
  prob = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=[tau] )
  msol2 =  solve( prob,  solver, callback=cb, saveat=dt )  
  plot( msol2, ; legend=false, xlim=(1999,2021), label="dde, with hsa, no fishing" )
end



