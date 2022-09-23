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


    
no_digits = 3  # time floating point rounding 

smallnumber = 1.0e-9 # floating point value sufficient to assume 0 valued
    


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

if aulab=="cfanorth" 
   ki = 1

elseif aulab=="cfasouth"
   ki = 2

elseif aulab=="cfa4x"
   ki = 3

end

kmu  = kmu = Kmu[ki] / mean(scale_factor)


# spin up time of ~ 1 cycle prior to start of dymamics and project 5 years into the future
tspan = (minimum(yrs) - 5.0, maximum(yrs) + nP + 1.1 )  

dt = 0.01 # time resolution of differ eq model solutions

survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)    # time of observations for survey
# survey_time = Y[:,:yrs]

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
hsa = Interpolations.scale(efc, yrs, 1:nS )

fish_time =  round.( round.( removals[:,:ts] ./ dt; digits=0 ) .* dt; digits=no_digits)    # time of observations for survey
# fish_time =  removals[:,:ts]     # time of observations for survey


removed = removals[:,Symbol("$aulab")]

function affect_fishing!(integrator)
  i = findall(t -> t == integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
end

# cb = CallbackSet( 
#   PresetTimeCallback( fish_time, affect_fishing! ), 
#   PositiveDomain()
# );

cb = PresetTimeCallback( fish_time, affect_fishing! ) 

# history function 0.5 default
# h(p,t) = ones( nS ) .* 0.5  #values of u before t0
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* kmu
 
tau = 1  # delay resolution
lags = [tau]

# allsavetimes = unique( vcat( survey_time, prediction_time  ) )

# stiff solvers: Rodas4()  ; Rosenbrock23()
# solver = MethodOfSteps(Rosenbrock23()) # slow  
# solver = MethodOfSteps(Rodas4())  
# other solvers: BS3() and Vern6() also RK4()
# solver = MethodOfSteps(Rodas5())  # safer 

# relative timings:
# solver = MethodOfSteps(Tsit5())  # 10 - 41.43 
# solver = MethodOfSteps(Rodas5())   # 20.94  - 71.73
# solver = MethodOfSteps(BS3())   # 56.1
# solver = MethodOfSteps(Rodas4()) #   24.86- 82.79
# solver = MethodOfSteps(Rosenbrock23()) #  71.48
# solver = MethodOfSteps(Vern6())  # 73.98
# solver = MethodOfSteps(RK4())   # 76.28
# solver = MethodOfSteps(TRBDF2())  # 92.16
# solver = MethodOfSteps(QNDF())  # 110.79
# solver = MethodOfSteps(Vern7())  #  111.7
# solver = MethodOfSteps(KenCarp4())  # 139.88


# solver = MethodOfSteps(Tsit5())  
solver = MethodOfSteps(Rodas5())  # safer 

Turing.setadbackend(:forwarddiff)  # only AD that works right now
Turing.setrdcache(true)

# these are dummy initial values .. just to get things started

u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* kmu
b=[1.0, 0.8]
K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu; 
d=[0.2, 0.3, 0.4, 0.5, 0.5, 0.5];
v=[0.8, 1.0, 1.0, 1.0];  
tau=1.0; 

p = ( b, K, d, v, tau, hsa )    # dummy values needs to start the turing initialization

if false
  ## test, can ignore
  prob = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=lags )
  msol2 =  solve( prob,  solver, callback=cb, saveat=dt )  
  plot( msol2, ; legend=false, xlim=(1999,2021), label="dde, with hsa, no fishing" )
end



