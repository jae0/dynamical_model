
using DifferentialEquations

if false
  # limit lower bound to positives
  function size_structured_dde!( du, u, h, p, t)
    uu = max.( 0.0, u )
    (b, K, d, v, tau, hsa)  = p
    tr21 = v[1] * max(0.0, h(p, t-1)[2])   # transition 2 -> 1   
    tr32 = v[2] * max(0.0, h(p, t-1)[3])   # transition 3 -> 2
    tr43 = v[3] * max(0.0, h(p, t-1)[4])   # transition 4 -> 3
    tr54 = v[4] * max(0.0, h(p, t-1)[5])   # transition 5 -> 4
    FP  = max(0.0, h(p, t-8)[6]   )      # no mature fem 8  yrs ago
    du[1] = tr21             - d[1] * uu[1] * ( uu[1] / hsa(t,1) )  # second order mortality       
    du[2] = tr32      - tr21 - d[2] * uu[2] * ( uu[2] / hsa(t,2) )  
    du[3] = tr43      - tr32 - d[3] * uu[3] * ( uu[3] / hsa(t,3) ) 
    du[4] = tr54      - tr43 - d[4] * uu[4] * ( uu[4] / hsa(t,4) ) 
    du[5] = b[1] * FP - tr54 - d[5] * uu[5] * ( uu[5] / hsa(t,5) )  
    du[6] = b[2] * FP        - d[6] * uu[6] * ( uu[6] / hsa(t,6) )   # fem mat simple logistic with lag tau and density dep on present numbers
  end

end


function size_structured_dde!( du, u, h, p, t)
    (b, K, d, v, tau, hsa)  = p
    # tr21 = v[1] * h(p, t-1)[2]    # transition 2 -> 1   
    # tr32 = v[2] * h(p, t-1)[3]    # transition 3 -> 2
    # tr43 = v[3] * h(p, t-1)[4]    # transition 4 -> 3
    # tr54 = v[4] * h(p, t-1)[5]    # transition 5 -> 4
    tr  = v .* h(p, t-1)[2:5]
    bu8 = b .* h(p, t-8)[6]      # no mature fem 8  yrs ago
    duh = d .* u ./ hsa(t,1:6)
    du[1] = tr[1]             - duh[1]     # note:       u = U /K; we are modifying K-> K*H; so d*u/H 
    du[2] = tr[2]     - tr[1] - duh[2]   
    du[3] = tr[3]     - tr[2] - duh[3]  
    du[4] = tr[4]     - tr[3] - duh[4]  
    du[5] = bu8[1]    - tr[4] - duh[5]   
    du[6] = bu8[2]            - duh[6]    # fem mat simple logistic with lag tau and density dep on present numbers
end

