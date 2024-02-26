# automatic load of things related to all projects go here

current_directory =  @__DIR__() 
print( "Current directory is: ", current_directory, "\n\n" )

pkgs_startup = [  
    "Revise", "Logging", "OhMyREPL",
    "Setfield", "Memoization",
    "DataFrames", "CSV", "RData",    
    "StatsBase", "Statistics", "Distributions", "Random", 
    "PlotThemes", "Colors", "ColorSchemes", 
    "Plots", "StatsPlots" 
]


# load directly can cause conflicts due to same function names 
pkgtoskipload = ["CairoMakie", "CairoMakie", "PlotlyJS",  "PlotlyBase",  "PlotlyKaleido" ]


if @isdefined pkgs 
    pkgs = unique!( [pkgs_startup; pkgs] )
else 
    pkgs = pkgs_startup
end


print( "Loading libraries:\n\n" )

# load libs and check settings
# pkgs are defined in snowcrab_startup.jl
using Pkg
for pk in pkgs; 
    if Base.find_package(pk) === nothing
        Pkg.add(pk)
    end
end   # Pkg.add( pkgs ) # add required packages

for pk in pkgs; 
    if !(Base.find_package(pk) === nothing)
        if !(pk in pkgtoskipload)
            @eval using $(Symbol(pk)); 
        end
    end
end



# functions required before loading environment

function discretize_decimal( x, delta=0.01 ) 
    num_digits = Int(ceil( log10(1.0 / delta)) )   # time floating point rounding
    out = round.( round.( x ./ delta; digits=0 ) .* delta; digits=num_digits)
    return out
end
  
colorscheme!("Monokai24bit") # for REPL

print( "\n\n", "Additional functions: \n\n" )

