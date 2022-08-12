 
 
# ------------------------------
# Delay SSA size structured

# https://github.com/palmtree2013/DelaySSAToolkit.jl
# https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/
# NOTE::: require 03.snowcrab_carstm.r to be completed 


# from https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/:

# "A major assumption behind the majority of stochastic models of biochemical kinetics is the memoryless hypothesis, i.e., the stochastic dynamics of the reactants is only influenced by the current state of the system, which implies that the waiting times for reaction events obey exponential distributions. Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics for such systems [1]. While this Markovian assumption considerably simplifies model analysis, it is dubious for modelling certain non-elementary reaction events that encapsulate multiple intermediate reaction steps [2]."

# For a few number of jumps, DelayRejection and DelayDirect will often perform better than other aggregators.

# For large numbers of jumps with sparse chain like structures and similar jump rates, for example continuous time random walks, DelayDirectCR and DelayMNRM often have the best performance.

 
pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions", "DynamicPPL",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  "Interpolations",
  "Plots", "StatsPlots", "MultivariateStats", "Graphviz_jll", "DelaySSAToolkit", "JumpProcesses"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end



rn = @reaction_network begin
  v1 * S2,             S2 --> S1
  v2 * S3,             S3 --> S2
  v3 * S4,             S4 --> S3
  v4 * S5,             S5 --> S4
  b1 * S6,             S6 --> S5
  b2 * S6,             S6 --> S6
end v1 v2 v3 v4 b1 b2 

species(rn)
parameters(rn)
Graph(rn)


tspan = (0.0, 100.0)
dt = 0.1

u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 100

b=[ 3.0, 1.24]
K=[100.0, 100.0, 100.0, 100.0, 100.0, 100.0];
d=[0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
v=[0.8, 0.8, 0.8, 0.8];  
p = ( 
  b1=b[1], b2=b[2],
  K1=K[1], K2=K[2], K3=K[3], K4=K[4], K5=K[5], K6=K[6], 
  d1=d[1], d2=d[2], d3=d[3], d4=d[4], d5=d[5], d6=d[6], 
  v1=v[1], v2=v[2], v3=v[3], v4=v[4] 
)

# instantaneous probability per time a jump occurs when
# the current state is u, current parameters are p, and the time is t
b1_rate(u,p,t) = p.b1
b1_affect!(integrator) = (integrator.u[1] += 1)
b1_jump = ConstantRateJump(b1_rate, b1_affect!)

d1_rate(u,p,t) = p.d1 * u[1]
d1_affect!(integrator) = (integrator.u[1] -= 1; integrator.u[2] += 1)
d1_jump = ConstantRateJump(d1_rate, d1_affect!)

h1_rate(u,p,t) = hsa(t,1) 
h1_affect!(integrator) = (integrator.u[1] -= 1);
h1_jump = VariableRateJump(h1_rate, h1_affect!)

# this is fora comre complex random rate (nonconstant rate) jump .. 
rng = JumpProcesses.DEFAULT_RNG
Random.TaskLocalRNG()
b1r_rate(u,p,t) = p.b1
# define the affect function via a closure

b1r_affect! = integrator -> let rng=rng
    # N(t) <-- N(t) + 1
    integrator.u[1] += 1
    # G(t) <-- G(t) + C_{N(t)}
    integrator.u[2] += rand(rng, (-1,1))
    nothing
end
b1r_jump = ConstantRateJump(b1r_rate, b1r_affect!)


dprob = DiscreteProblem( u0, tspan, p )
jprob = JumpProblem(dprob, Direct(), b1_jump, h1_jump, d1_jump; save_positions = (false, false) )
sol = solve(jprob, SSAStepper(); saveat=dt )
plot(sol, label="N(t)", xlabel="t", legend=:bottomright)

 

delay_trigger_affect1! = function (integrator, rng)
  append!(integrator.de_chan[1], 1)
end
delay_trigger_affect2! = function (integrator, rng)
  append!(integrator.de_chan[2], 1)
end
delay_trigger_affect3! = function (integrator, rng)
  append!(integrator.de_chan[3], 1)
end
delay_trigger_affect4! = function (integrator, rng)
  append!(integrator.de_chan[4], 1)
end
delay_trigger_affect5! = function (integrator, rng)
  append!(integrator.de_chan[5], 8)
end
delay_trigger_affect6! = function (integrator, rng)
  append!(integrator.de_chan[6], 8)
end
 
delay_trigger = Dict(
  7=>delay_trigger_affect1!,
  8=>delay_trigger_affect2!,
  9=>delay_trigger_affect3!,
  10=>delay_trigger_affect4!,
  11=>delay_trigger_affect5!,
  12=>delay_trigger_affect6!
)

delay_interrupt = Dict()
delay_complete = Dict(
  1=>[1=>1, 2=>-1],
  2=>[2=>1, 3=>-1],
  3=>[3=>1, 4=>-1],
  4=>[4=>1, 5=>-1],
  5=>[5=>1],
  6=>[6=>1]
)

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

de_chan0 = [[]]

djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0 )
  
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, 
  saveat=dt, save_positions=(false,false))

sol = solve(djprob, SSAStepper())


function size_structured!( du, u, h, p, t)
  # S1, S2, S3, S4, S5, S6 = u
  b, K, d, v, tau, hsa  = p
  tr21 = v[1] * h(p, t-1)[2]   # transition 2 -> 1   
  tr32 = v[2] * h(p, t-1)[3]   # transitiom 3 -> 2
  tr43 = v[3] * h(p, t-1)[4]   # transitiom 4 -> 3
  tr54 = v[4] * h(p, t-1)[5]   # transitiom 5 -> 4
  f8  = h(p, t-8)[6]           # no fem 8 yrs ago
  du[1] = tr21             - (d[1] * u[1]) * (u[1]/ K[1]) * hsa(t,1)       
  du[2] = tr32      - tr21 - (d[2] * u[2]) * (u[2]/ K[2]) * hsa(t,2) 
  du[3] = tr43      - tr32 - (d[3] * u[3]) * (u[3]/ K[3]) * hsa(t,3)
  du[4] = tr54      - tr43 - (d[4] * u[4]) * (u[4]/ K[4]) * hsa(t,4)
  du[5] = b[1] * f8 - tr54 - (d[5] * u[5]) * (u[5]/ K[5]) * hsa(t,5) 
  du[6] = b[2] * f8        - (d[6] * u[6]) * (u[6]/ K[6]) * hsa(t,6)  # fem mat simple logistic with lag tau and density dep on present numbers
end



  
