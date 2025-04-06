using Random
using DifferentialEquations
using Plots
using Distributions
using RotaryInvertedPendulum
using LinearAlgebra

#create model
include("parameters.jl")
(M,N,M_f,N_f), sys_cont, Total_energy, T_f,V_f = generate_dynamics(dyn_params,x_equil)

Damping_dist = truncated(Normal(Damping, 0.0001), 0.5*Damping, 2*Damping)
trajectories = 100

# Generate control input from polynomial coefficients
function polynomial_control(coeffs, t,T)
    degree = length(coeffs) - 1
    return sum(coeffs[i+1] * (t/T)^i for i in 0:degree)
end

n_coeffs=10

# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl_test!(du, u, p, t)
    θ=collect(u[1:2])
    θd=collect(u[3:4])
    
    Damping_ratio=p[1]
    coeffs=p[2:n_coeffs+1]

    u_f(t)=polynomial_control(coeffs, t,1)

    D=[0 0;0 Damping_ratio]
    #D=[0 0;0 0]
    Damping_force=D*θd
  
    M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ...),N_f(u...),Damping_force)
    
     # Fill in-place the du vector (du = dx/dt = [θd; θ̈])
     du[1:2] .= θd
     du[3:4] .= M_a \ (B_a * u_f(t) - N_a)
    #return vec([θd;inv(M_a)*(B_a*u_f(t)-N_a)])
end

#simulate!
q0=[0.1;pi-0.2;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/10;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, 5.0)
ps=[rand(Random.GLOBAL_RNG,Damping_dist);0.01*ones(n_coeffs)]
prob = ODEProblem(dynamics_acc_ctrl_test!, q0, tspan, ps)
sol = solve(prob, Tsit5())
#Simulate & animate!
#tvec,q_sol,qd_sol=pend_sim(prob)#,reltol=1e-10,abstol=1e-10)
#plot(tvec,q_sol)

# Define ensemble with randomized damping
ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> begin
    new_p = [rand(Damping_dist);ps[2:end]]
    remake(prob, p = new_p)
end)

# Solve
ensemblesol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 100)

using Optimization, OptimizationNLopt, OptimizationMOI, SciMLExpectations

gd = GenericDistribution(Damping_dist)

#h(x, u, p) = u, [p[1];x[1]]
obs(sol, p) = abs2(sol[2, end] - pi) + abs2(sol[3, end]) + abs2(sol[4, end])
#h(x,u, p) =u, p
h(x,u, p) = u, p

function 𝔼_loss(coeffs, pars)
    full_p = [Damping; coeffs]
    prob = ODEProblem(dynamics_acc_ctrl_test!, q0, tspan, full_p)
    sm = SystemMap(prob, Tsit5())
    exprob = ExpectationProblem(sm, obs, h, gd)
    sol = solve(exprob, Koopman(), ireltol=1e-5)
    return sol.u  # scalar loss
end

coeffs_init = 0.3*ones(n_coeffs)
# p = [0.3] is the fixed value passed to the function; just a placeholder
opt_f = OptimizationFunction(𝔼_loss, Optimization.AutoForwardDiff())
opt_ini = 0.3*ones(n_coeffs)
opt_lb = -0.1*ones(n_coeffs)
opt_ub = 100*ones(n_coeffs)
opt_prob = OptimizationProblem(opt_f, coeffs_init, p=[Damping],
                               lb = -1000 * ones(n_coeffs),
                               ub = 1000 * ones(n_coeffs))
                               optimizer = OptimizationMOI.MOI.OptimizerWithAttributes(
                                NLopt.Optimizer, "algorithm" => :LD_MMA
                            )
            

optimizer = OptimizationMOI.MOI.OptimizerWithAttributes(NLopt.Optimizer,
    "algorithm" => :LD_MMA)
opt_sol = solve(opt_prob, optimizer)
minx = opt_sol.u


# Define ensemble with randomized damping
ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> begin
    new_p = [rand(Damping_dist);minx]
    remake(prob, p = new_p)
end)

ensemblesol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 100)
