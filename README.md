# NOTE: This project is in maintenance mode. Further development (backwards incompatible) of more general models is occuring on: 

[https://bitbucket.org/autocatalysis/model_fishery](https://bitbucket.org/autocatalysis/model_fishery)

This github instance is here to serve as a reference for the online manuscript
on bioRxiv (see below).

# Dynamical models for fishery assessment

They are made possible by some excellent work:

- [Julia](https://julialang.org/)
- [SciML](https://sciml.ai/) 
- [Turing](https://turing.ml/) 
- [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
  
- Other important libraries used can be found in the set up environment (*.environment.jl) for each project.


---

## Current results:

- [Results based upon this developmental branch.](snowcrab/snowcrab_results_current.md)

- They are generated by the [main calling program.](snowcrab/04.snowcrab_fishery_model.md) 


---

## Models currently implemented:

- Discrete logistic (Euler form) 

- Discrete logistic map 

- Six sex/size-stage continuous delay differential equation model  



## Example data are in RData form: 

- [biodyn_number_size_struct.RData](snowcrab/data/biodyn_number_size_struct.RData) 
        
    - time series by large areal units, number, size, habitat, etc.

    - used by stage structured models

- [biodyn_biomass.RData](snowcrab/data/biodyn_biomass.RData) 
        
    - time series by large areal units in biomass 

    - used by logistic models 


---

## Documentation:

[An example using this approach is found here](https://www.biorxiv.org/content/10.1101/2023.02.13.528296v3) for **snowcrab**. It is now grafted onto bio.snowcrab for operational work (as of Jan 2023) 

[Current documentation is here.](snowcrab/ms_choi_v3.pdf)  
