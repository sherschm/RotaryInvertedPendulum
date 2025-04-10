using Random
using DifferentialEquations
using Plots
using Distributions
using RotaryInvertedPendulum
using LinearAlgebra
using DelimitedFiles

#create model
include("parameters.jl")
(M,N,M_f,N_f), sys_cont, Total_energy, T_f,V_f = generate_dynamics(dyn_params,x_equil)

Damping_dist = truncated(Normal(Damping, 0.00005), 0.5*Damping, 2*Damping)
trajectories = 70

# Generate control input from polynomial coefficients
function polynomial_control(coeffs, t,T)
    degree = length(coeffs) - 1
    return sum(coeffs[i+1] * (t/T)^i for i in 0:degree)
end

function chebyshev_control(t, coeffs, T_final)
    # Time normalization: [0, T_final] → [-1, 1]
    τ = 2t / T_final - 1

    # Evaluate Chebyshev series at normalized time
    N = length(coeffs)
    T = [cos(i * acos(τ)) for i in 0:N-1]  # Chebyshev basis at τ

    return dot(coeffs, T)
end

function bernstein_control(t, coeffs, T_final)
    n = length(coeffs) - 1  # degree of the polynomial
    τ = t / T_final         # normalize time to [0, 1]

    u = 0.0
    for i in 0:n
        B = binomial(n, i) * (1 - τ)^(n - i) * τ^i
        u += coeffs[i + 1] * B
    end

    return u
end

Tf=6
n_coeffs=50

# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl_test!(du, u, p, t)
    θ=collect(u[1:2])
    θd=collect(u[3:4])
    
    Damping_ratio=p[1]
    coeffs=p[2:n_coeffs+1]

    #u_f(t)=polynomial_control(coeffs, t,Tf)
    u_f(t)=chebyshev_control(t,coeffs, Tf)
   # u_f(t)=bernstein_control(t,coeffs, Tf)

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
q0=[0.0;0.0;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/10;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, Tf)
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
#obs(sol, p) = abs2(sol[2, end] - pi) + abs2(sol[3, end]) + abs2(sol[4, end])

λ = 0.001  # regularization strength
#=obs(sol, p) = begin
    coeffs = p[2:end]  # assuming p[1] = damping, p[2:end] = Chebyshev coeffs
    reg_term = λ * sum(abs2, coeffs)
    terminal_cost = 10*abs2(sol[2, end] - π) + 0.1*abs2(sol[4, end]) + 0.1*abs2(sol[3, end]) + 0.1*abs2(sol[1, end])# + abs2(sol[3, end]) + abs2(sol[4, end])

    return terminal_cost + reg_term
end# =#

obs(sol, p) = begin
    coeffs = p[2:end]
    terminal_cost = 100 * abs2(sol[2, end] - (π+0.1)) + abs2(sol[3, end]) + abs2(sol[1, end])+abs2(sol[4, end])
    
    reg_term = λ * sum(abs2, coeffs)
    dt = sol.t[end] / (length(sol.t) - 1)
    # Running cost to penalize being far from π throughout the trajectory
   #run_cost = sum(abs2, + sol[3, :]) * dt   # integrate over time
   # return terminal_cost +  0.1*run_cost + reg_term
    return terminal_cost + reg_term
end
#=
obs(sol, p) = begin
    coeffs = p[2:end]
    terminal_cost = 10 * sum(abs2,(sol[2, end-5:end] .- π)) #+abs2(sol[4, end])
    
    reg_term = λ * sum(abs2, coeffs)
    return terminal_cost  + reg_term
end
=#
#h(x,u, p) =u, p
#h(x,u, p) = u, [p[1];x[1]]

h(x,u, p) = u, p
prob_template = ODEProblem(dynamics_acc_ctrl_test!, q0, tspan, ps)

function 𝔼_loss(coeffs,pars)

    full_p=[pars[1];coeffs]
   # Damping = full_p[1]
    #coeffs = full_p[2:end]

    prob = remake(prob_template, p = full_p)
    sm = SystemMap(prob, Tsit5())
    exprob = ExpectationProblem(sm, obs, h, gd)
    sol = solve(exprob, Koopman(), ireltol=1e-5)
    return sol.u  # scalar loss
end

coeffs_init = zeros(n_coeffs)
opt_f = OptimizationFunction(𝔼_loss, Optimization.AutoForwardDiff())

full_p_init = vcat(Damping,coeffs_init)
opt_prob = OptimizationProblem(opt_f, coeffs_init, Damping)
                               #lb = -1000 * ones(n_coeffs),
                               #ub = 1000 * ones(n_coeffs)
                               
        
optimizer = OptimizationMOI.MOI.OptimizerWithAttributes(NLopt.Optimizer,
    "algorithm" => :LD_MMA)
opt_sol = solve(opt_prob, optimizer)
optimal_coefficients = opt_sol.u
println(optimal_coefficients)

# Define ensemble with randomized damping
ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> begin
    new_p = [rand(Damping_dist);optimal_coefficients]
    remake(prob, p = new_p)
end)

ensemblesol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 100)

prob = ODEProblem(dynamics_acc_ctrl_test!, q0, tspan, [rand(Damping_dist);optimal_coefficients])
#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(prob)#,reltol=1e-10,abstol=1e-10)
plot(tvec,q_sol)
plot_params=(l1,l2)
rot_pendulum_animator(q_sol,tvec,plot_params;name="swing_up_robust")

tvec_out=0.005:0.005:Tf
acc_cmd_out=zeros(length(tvec_out))

for i in 1:length(tvec_out)

    acc_cmd_out[i]=chebyshev_control(tvec_out[i], optimal_coefficients, Tf)
end

writedlm("data/swing_up/swingup_acc_robust_cmd.csv", [Float32.(tvec_out) Float32.(acc_cmd_out/encoder_steps_per_rad)],",")
