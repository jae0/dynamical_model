
using DifferentialEquations

if false
function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  
    u = max.(u, 1.0)
    (b, K, d, v, tau, hsa)  = p
    u1 = max.(h(p, t-1.0), 1.0)   # no in previous years
    u8 = max.(h(p, t-8.0), 1.0)  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6) 
    # NOTE:: breaking it down like this actually helps to make it faster  ;)
    du[1] = v[1] * u1[2]                 - d[1] * u[1] * u[1]  / K[1] / vh[1]       # note:       
    du[2] = v[2] * u1[3]  - v[1] * u1[2] - d[2] * u[2] * u[2]  / K[2] / vh[2]     
    du[3] = v[3] * u1[4]  - v[2] * u1[3] - d[3] * u[3] * u[3]  / K[3] / vh[3]    
    du[4] = v[4] * u1[5]  - v[3] * u1[4] - d[4] * u[4] * u[4]  / K[4] / vh[4]    
    du[5] = b[1] * u8[6]  - v[4] * u1[5] - d[5] * u[5] * u[5]  / K[5] / vh[5]     
    du[6] = b[2] * u8[6]                 - d[6] * u[6] * u[6]  / K[6] / vh[6]      # fem mat simple logistic with lag tau and density dep on present numbers
end


function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  
    u = max.(u, 1.0)
    (b, K, d, v, tau, hsa)  = p
    u1 = max.(h(p, t-1.0), 1.0)   # no in previous years
    u8 = max.(h(p, t-8.0), 1.0)  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6) 
    # NOTE:: breaking it down like this actually helps to make it faster  ;)
    du[1] = v[1] * u1[2]                 - d[1] * u[1] /  vh[1]       # note:       
    du[2] = v[2] * u1[3]  - v[1] * u1[2] - d[2] * u[2] /  vh[2]     
    du[3] = v[3] * u1[4]  - v[2] * u1[3] - d[3] * u[3] /  vh[3]    
    du[4] = v[4] * u1[5]  - v[3] * u1[4] - d[4] * u[4] /  vh[4]    
    du[5] = b[1] * u8[6]  - v[4] * u1[5] - d[5] * u[5] /  vh[5]     
    du[6] = b[2] * u8[6]                 - d[6] * u[6] /  vh[6]      # fem mat simple logistic with lag tau and density dep on present numbers
end

function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  
    u = max.(u, 1.0)
    (b, K, d, v, tau, hsa)  = p
    u1 = max.(h(p, t-1.0), 1.0)   # no in previous years
    u8 = max.(h(p, t-8.0), 1.0)  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6) 
    # NOTE:: breaking it down like this actually helps to make it faster  ;)
    du[1] = v[1] * u1[2]                 - d[1] * u[1]           # note:       
    du[2] = v[2] * u1[3]  - v[1] * u1[2] - d[2] * u[2]         
    du[3] = v[3] * u1[4]  - v[2] * u1[3] - d[3] * u[3]        
    du[4] = v[4] * u1[5]  - v[3] * u1[4] - d[4] * u[4]        
    du[5] = b[1] * u8[6]  - v[4] * u1[5] - d[5] * u[5]         
    du[6] = b[2] * u8[6]                 - d[6] * u[6]          # fem mat simple logistic with lag tau and density dep on present numbers
end


function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  
    u = max.(u, 0.0)
    (b, K, d, dh, v, tau, hsa)  = p
    u1 = max.(h(p, t-1.0), 0.0)   # no in previous years
    u8 = max.(h(p, t-8.0), 0.0)  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6)  
    # NOTE:: breaking it down like this seems to actually helps to make it faster  ;)
    du[1] = v[1] * u1[2]                 - d[1] * u[1] -  dh[1] * u[1] * u[1] / K[1] / vh[1]       # note:       
    du[2] = v[2] * u1[3]  - v[1] * u1[2] - d[2] * u[2] -  dh[2] * u[2] * u[2] / K[2] / vh[2]     
    du[3] = v[3] * u1[4]  - v[2] * u1[3] - d[3] * u[3] -  dh[3] * u[3] * u[3] / K[3] / vh[3]    
    du[4] = v[4] * u1[5]  - v[3] * u1[4] - d[4] * u[4] -  dh[4] * u[4] * u[4] / K[4] / vh[4]    
    du[5] = b[1] * u8[6]  - v[4] * u1[5] - d[5] * u[5] -  dh[5] * u[5] * u[5] / K[5] / vh[5]     
    du[6] = b[2] * u8[6]                 - d[6] * u[6] -  dh[6] * u[6] * u[6] / K[6] / vh[6]      # fem mat simple logistic with lag tau and density dep on present numbers
end

end


function size_structured_dde!( du, u, h, p, t)
    # here u, du are actual numbers .. not normalized by K due to use of callbacks  
    u = max.(u, 0.0)
    (b, K, d, dh, v, tau, hsa)  = p
    u1 = max.(h(p, t-1.0), 0.0)   # no in previous years
    u8 = max.(h(p, t-8.0), 0.0)  # no mature fem 8  yrs ago
    vh = hsa(t, 1:6)  
    # NOTE:: breaking it down like this seems to actually helps to make it faster  ;)
    du[1] = v[1] * u1[2]                 - d[1] * u[1] * u[1] / K[1] -  dh[1] * u[1] * u[1] / K[1] / vh[1]       # note:       
    du[2] = v[2] * u1[3]  - v[1] * u1[2] - d[2] * u[2] * u[2] / K[2] -  dh[2] * u[2] * u[2] / K[2] / vh[2]     
    du[3] = v[3] * u1[4]  - v[2] * u1[3] - d[3] * u[3] * u[3] / K[3] -  dh[3] * u[3] * u[3] / K[3] / vh[3]    
    du[4] = v[4] * u1[5]  - v[3] * u1[4] - d[4] * u[4] * u[4] / K[4] -  dh[4] * u[4] * u[4] / K[4] / vh[4]    
    du[5] = b[1] * u8[6]  - v[4] * u1[5] - d[5] * u[5] * u[5] / K[5] -  dh[5] * u[5] * u[5] / K[5] / vh[5]     
    du[6] = b[2] * u8[6]                 - d[6] * u[6] * u[6] / K[6] -  dh[6] * u[6] * u[6] / K[6] / vh[6]      # fem mat simple logistic with lag tau and density dep on present numbers
end
