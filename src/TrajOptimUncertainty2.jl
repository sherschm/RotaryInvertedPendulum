using Random, Distributions, Optim, Plots, DifferentialEquations

# Define the inverted pendulum dynamics as an ODE
function pendulum_dynamics!(du, u, p, t)
    x, θ, x_dot, θ_dot = u
    coeffs, T = p
    control = sum(coeffs[i+1] * (t/T)^i for i in 0:(length(coeffs)-1))
    
    g = 9.81  # gravity (m/s^2)
    m = 1.0   # pendulum mass (kg)
    l = 1.0   # pendulum length (m)
    b = 0.1   # damping coefficient
    I = m * l^2  # moment of inertia
    
    θ_ddot = (m * g * l * sin(θ) - b * θ_dot + control) / I
    x_ddot = control / m
    
    du[1] = x_dot
    du[2] = θ_dot
    du[3] = x_ddot
    du[4] = θ_ddot
end

# Define the cost function
function cost_function(coeffs)
    dt = 0.01  # time step
    T = 5.0    # total time (longer for swing-up)
    tspan = (0.0, T)
    u0 = [0.0, 0.0, 0.0, 0.0]  # initial state (x, θ, x_dot, θ_dot)
    p = (coeffs, T)
    prob = ODEProblem(pendulum_dynamics!, u0, tspan, p)
    sol = solve(prob, Tsit5(), dt=dt)
    
    #cost = sum((θ - π)^2 + 0.01 * u^2 for (θ, u) in zip(sol[2, :], sol[3, :]))  # Penalize deviation and control effort
    cost = sum((θ - π)^2 + 0.01 * u^2 for (θ, u) in zip(sol[2, :], sol[3, :]))  # Penalize deviation and control effort
    return cost
end

# Optimization setup
function optimize_trajectory()
    degree = 10  # Polynomial degree for control function
    initial_coeffs = ones(degree + 1)
    
    result = optimize(cost_function, initial_coeffs, LBFGS())
    return result.minimizer
end

# Run the optimization and animate results
optimal_coeffs = optimize_trajectory()
println("Optimal coefficients:", optimal_coeffs)
#animate_pendulum(optimal_coeffs)

# Simulate and animate the optimized response
#function animate_pendulum(coeffs)
dt = 0.05
T =10.0
tspan = (0.0, T)
u0 = [0.0, 0.0, 0.0, 0.0]
p = (optimal_coeffs, T)
prob = ODEProblem(pendulum_dynamics!, u0, tspan, p)
sol = solve(prob, Tsit5(), dt=dt)

anim = Animation()
for i in 1:length(sol.t)
    θ = sol[2, i]
    x = sin(θ)  # pendulum end x
    y = -cos(θ) # pendulum end y
    
    plt = plot([-x, x], [0, y], xlims=(-1.2, 1.2), ylims=(-1.2, 1.2), linewidth=3, label=false,aspect_ratio=:equal)
    frame(anim, plt)
end
gif(anim, "pendulum_animation.gif", fps=20)
#end