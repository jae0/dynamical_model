# Dynamical models for fishery assessment made possible by [Julia](https://julialang.org/), [SciML](https://sciml.ai/) and [Turing](https://turing.ml/). 


---

[Current results based upon this developmental branch are found here.](snowcrab/snowcrab_results_current.md)

---

An example using this approach is found in (https://www.biorxiv.org/content/10.1101/2023.02.13.528296v3) for **snowcrab** -- now grafted onto bio.snowcrab for operational work (Jan 2023) 

---

This instance is a developmental branch. Current models:

- discrete logistic (Euler) (logistic_discrete\*)

- discrete logistic map (logistic_discrete\*)

- six size-stage continuous delay difference (size_structured_dde\*)

- [documentation](https://www.biorxiv.org/content/10.1101/2023.02.13.528296v3)  



- Example data in RData form: 

    - biodyn_number_size_struct.RData (time series by large areal units, number, size, habitat, etc)

    - biodyn_biomass.RData (time series by large areal units in biomass)
    
- [Main program is here.](https://github.com/jae0/dynamical_model/blob/master/snowcrab/README.md), adjust file paths as required

- Julia startup files define environment, adjust as required.

