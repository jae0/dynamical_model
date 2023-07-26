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

if @isdefined pkgs 
    pkgs = unique!( [pkgs_startup; pkgs] )
else 
    pkgs = pkgs_startup
end
    
print( "Loading libraries:\n\n" )
for pk in pkgs; 
    print(pk, "\n")
    @eval using $(Symbol(pk)); 
end    


# functions required before loading environment

function discretize_decimal( x, delta=0.01 ) 
    num_digits = Int(ceil( log10(1.0 / delta)) )   # time floating point rounding
    out = round.( round.( x ./ delta; digits=0 ) .* delta; digits=num_digits)
    return out
end
  
colorscheme!("Monokai24bit") # for REPL

print( "\n\n", "Additional functions: \n\n" )

