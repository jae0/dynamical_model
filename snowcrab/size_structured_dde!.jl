
using DifferentialEquations
 
function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  

    (b, K, d, v, tau, hsa)  = p
    u1 = h(p, t-1.0)    # no in previous years
    f8 = h(p, t-8.0)[6]  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6)  
    
    br =  f8 .* b    
    tr =  v .* u1[2:5]
    dr =  d .* u .* u ./K ./vh 
    
    du[1] = tr[1]            - dr[1]       # note:       
    du[2] = tr[2]   - tr[1]  - dr[2]     
    du[3] = tr[3]   - tr[2]  - dr[3]    
    du[4] = tr[4]   - tr[3]  - dr[4]    
    du[5] = br[1]   - tr[4]  - dr[5]     
    du[6] = br[2]            - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers
    
    # i = findall( x -> x < 0, u .+ du ) 
    # du[i] .+= u[i] .+ smallnumber  # force u to be postive
    
end



