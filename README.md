# dynamical_model

Dynamical models for fishery assessment in Julia
(https://www.biorxiv.org/content/10.1101/2023.02.13.528296v1)

**snowcrab** (dynamical only): -- now grafted onto bio.snowcrab for operational work (Jan 2023) 

.. this is here as a seprate form for exploratory work (to add spatiotemporal struture)

- discrete logistic (Euler) (logistic_discrete\*)
- discrete logistic map (logistic_discrete\*)
- six size-stage continuous delay difference (size_structured_dde\*)

- documentation in tex and pdf form
    - ms_choi_v2.pdf

- example data in RData form: 
    - biodyn_number_size_struct.RData (time series by large areal units, number, s
ize, habitat, etc)
    - biodyn_biomass.RData (time series by large areal units in biomass)
    
- main program: 04.snowcrab_fishery_model.jl, adjust file paths as required

- Julia startup files define environment, adjust as required.



