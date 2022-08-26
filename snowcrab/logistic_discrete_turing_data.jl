
# logistic_discrete_turing_data

# perpare dat for dde run of fishery model
using DifferentialEquations, Interpolations, RData
 
fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_biomass.RData"
o = load( fndat, convert=true)
Y = o["Y"]
Ksd = o["Ksd"]
Kmu = o["Kmu"]
removals = o["L"] 

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
        # fishery_model_data_inputs( year.assessment=2021, type="biomass_dynamics",  for_julia=TRUE   )
        
        # # then back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
        # @rget Y  
        # @rget Kmu 
        # @rget removals 
        # @rget ty
 
    end


# "survey index"
S = Y[:,Symbol("$aulab"  )]

nT = length(S)
# nP = 1 # no year to project into future

# convert to (biomass kt to number) 
kmu = Kmu[au]  
ksd = Ksd[au]

survey_time = Y[:,:yrs]   # time of observations for survey
prediction_time = floor.(vcat( survey_time, collect(1:nP) .+ maximum(survey_time) ) ) .+ round( round( 9.0/12.0 /dt; digits=0 ) *dt; digits=no_digits)
fish_time = round.( round.( removals[:,:ts] ./ dt; digits=0 ) .* dt; digits=no_digits)   # time of observations for survey

removed = removals[:,Symbol("$aulab")]
 
 
