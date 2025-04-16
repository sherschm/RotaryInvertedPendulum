using Random
using DifferentialEquations
using Plots
using Distributions
using RotaryInvertedPendulum
using LinearAlgebra
using DelimitedFiles
using Optimization, OptimizationNLopt, OptimizationMOI, SciMLExpectations

#create model
include("parameters.jl")
(M,N,M_f,N_f), sys_cont, Total_energy, T_f,V_f = generate_dynamics(dyn_params,x_equil)

Damping_sd=0.00005
Damping_dist = truncated(Normal(Damping, Damping_sd), 0.5*Damping, 2*Damping)
#R_
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

function fourier_control(t::Real, coeffs::Vector{<:Real}, T::Real)
    N = div(length(coeffs)-1, 2)  # Number of harmonics
    a0 = coeffs[1]
    u = a0 / 2

    for n in 1:N
        an = coeffs[2n]
        bn = coeffs[2n + 1]
        u += an * cos(2π * n * t / T) + bn * sin(2π * n * t / T)
    end

    return u
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

control_parametrisation="fourier"
Tf=4
n_coeffs=100
#=
Tf=6
n_coeffs=50=#
if control_parametrisation == "fourier"
    control(t, coeffs, T_final)= fourier_control(t, coeffs, T_final)
elseif control_parametrisation == "chebyshev"
    control(t, coeffs, T_final)= chebyshev_control(t, coeffs, T_final)
elseif control_parametrisation == "bernstein"
    control(t, coeffs, T_final)= bernstein_control(t, coeffs, T_final)
end

# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl_test!(du, u, p, t)
    θ=collect(u[1:2])
    θd=collect(u[3:4])
    
    Damping_ratio=p[1]
    coeffs=p[2:n_coeffs+1]

    u_f(t)=control(t,coeffs, Tf)

    D=[0 0;0 Damping_ratio]
    #D=[0 0;0 0]
    Damping_force=D*θd
  
    M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ...),N_f(u...),Damping_force)
    
    du[1:2] .= θd
    du[3:4] .= M_a \ (B_a * u_f(t) - N_a)
end

#simulate!
q0=[0.0;0.0;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/10;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, Tf)
ps=[rand(Random.GLOBAL_RNG,Damping_dist);0.01*ones(n_coeffs)]
prob = ODEProblem(dynamics_acc_ctrl_test!, q0, tspan, ps)
sol = solve(prob, Tsit5())

# Define ensemble with randomized damping
ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> begin
    new_p = [rand(Damping_dist);ps[2:end]]
    remake(prob, p = new_p)
end)

# Solve 
ensemblesol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 100)

gd = GenericDistribution(Damping_dist)

λ = 0.001  # regularization strength
obs(sol, p) = begin
    coeffs = p[2:end]
    terminal_cost = 1000 * abs2(sol[2, end] - (π)) + abs2(sol[3, end]) + abs2(sol[1, end]) + 10*abs2(sol[4, end])+control(Tf, coeffs, Tf)
    #terminal_cost = 100 * abs2(sol[2, end-0.3] - (π+0.1)) + abs2(sol[3, end-0.3]) + abs2(sol[1, end-0.3])+abs2(sol[4, end-0.3])
    reg_term = λ * sum(abs2, coeffs)
    dt = sol.t[end] / (length(sol.t) - 1)
    
    # Running cost to penalize being far from π throughout the trajectory
   #run_cost = sum(abs2, + sol[3, :]) * dt   # integrate over time
    #return terminal_cost +  0.1*run_cost + reg_term
    return terminal_cost + reg_term
end

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
    acc_cmd_out[i]=control(tvec_out[i], optimal_coefficients, Tf)
end

data_dir="data/swing_up/robust/"
writedlm(data_dir*"swingup_acc_robust_cmd.csv", [Float32.(tvec_out) Float32.(acc_cmd_out/encoder_steps_per_rad)],",")

p1=plot(ensemblesol,(vars=1),xlabel="Time (s)",ylabel="Theta 1 (rad)",dpi=300)
p2=plot(ensemblesol,(vars=2),xlabel="Time (s)",ylabel="Theta 2 (rad)",dpi=300)
p3=plot(ensemblesol,(vars=3),xlabel="Time (s)",ylabel="Theta 1 velocity (rad/s)",dpi=300)
p4=plot(ensemblesol,(vars=4),xlabel="Time (s)",ylabel="Theta 2 velocity (rad/s)",dpi=300)

plot(p1,p2,p3,p4,layout=(2,2),size=(700,500),dpi=300)
savefig(data_dir*"trajectory_reponse")

plot(tvec_out,acc_cmd_out,xlabel="Time (s)", ylabel= "Acceleration (rad/s^2)")
savefig(data_dir*"ctrl_acceleration")

function write_variable_descriptions(filename::String)
    open(filename, "w") do io
        # Write header
        println(io, "Variable Descriptions")
        println(io, "======================\n")

        println(io, "Joint 2 damping mean = $Damping, standard devation = $Damping_sd")
        println(io, "Trajectory time = $Tf s.")

        println(io, "Control function = "*control_parametrisation)
        println(io, "Number of control coefficients = $n_coeffs")
    end

    println("Variable descriptions written to $filename")
end

write_variable_descriptions(data_dir*"trajectory_parameters.txt")