
using DifferentialEquations



function size_structured!( du, u, h, p, t)
  uu = max.( u, eps )
  (b, K, d, v, tau, hsa)  = p
  tr21 = v[1] * max(eps, h(p, t-1)[2])   # transition 2 -> 1   
  tr32 = v[2] * max(eps, h(p, t-1)[3])   # transitiom 3 -> 2
  tr43 = v[3] * max(eps, h(p, t-1)[4])   # transitiom 4 -> 3
  tr54 = v[4] * max(eps, h(p, t-1)[5])   # transitiom 5 -> 4
  FP  = max(eps, h(p, t-8)[6]   )      # no mature fem 8  yrs ago
  du[1] = tr21             - d[1] * uu[1] * (uu[1] / (K[1]*hsa(t,1)) )  # second order mortality       
  du[2] = tr32      - tr21 - d[2] * uu[2] * (uu[2] / (K[2]*hsa(t,2)) )  
  du[3] = tr43      - tr32 - d[3] * uu[3] * (uu[3] / (K[3]*hsa(t,3)) ) 
  du[4] = tr54      - tr43 - d[4] * uu[4] * (uu[4] / (K[4]*hsa(t,4)) ) 
  du[5] = b[1] * FP - tr54 - d[5] * uu[5] * (uu[5] / (K[5]*hsa(t,5)) )  
  du[6] = b[2] * FP        - d[6] * uu[6] * (uu[6] / (K[6]*hsa(t,6)) )   # fem mat simple logistic with lag tau and density dep on present numbers
end


function size_structured_nomax!( du, u, h, p, t)
  (b, K, d, v, tau, hsa)  = p
  uu = 
  tr21 = v[1] * h(p, t-1)[2]    # transition 2 -> 1   
  tr32 = v[2] * h(p, t-1)[3]    # transitiom 3 -> 2
  tr43 = v[3] * h(p, t-1)[4]    # transitiom 4 -> 3
  tr54 = v[4] * h(p, t-1)[5]    # transitiom 5 -> 4
  FP  = h(p, t-8)[6]       # no mature fem 8  yrs ago
  du[1] = tr21             - d[1] * u[1] * (u[1] / (K[1]*hsa(t,1)) )  # second order mortality       
  du[2] = tr32      - tr21 - d[2] * u[2] * (u[2] / (K[2]*hsa(t,2)) )  
  du[3] = tr43      - tr32 - d[3] * u[3] * (u[3] / (K[3]*hsa(t,3)) ) 
  du[4] = tr54      - tr43 - d[4] * u[4] * (u[4] / (K[4]*hsa(t,4)) ) 
  du[5] = b[1] * FP - tr54 - d[5] * u[5] * (u[5] / (K[5]*hsa(t,5)) )  
  du[6] = b[2] * FP        - d[6] * u[6] * (u[6] / (K[6]*hsa(t,6)) )   # fem mat simple logistic with lag tau and density dep on present numbers
end
