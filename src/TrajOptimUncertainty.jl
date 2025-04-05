using Random
using DifferentialEquations
using Plots

# Define the implicit inverted pendulum dynamics
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

# Generate control input from polynomial coefficients
function polynomial_control(coeffs, t, T)
    degree = length(coeffs) - 1
    return sum(coeffs[i+1] * (t/T)^i for i in 0:degree)
end

# Cost function (to minimize time to reach an upright position)
function cost_function(solution)
    final_time = solution.t[end]
    final_angle = solution.u[end][1]
    final_velocity = solution.u[end][2]
    
    # Penalize if we don't end up upright (θ ~ 0) or not at rest (θ_dot ~ 0)
    #cost = final_time + 10 * abs(final_angle) + 1000 * abs(final_velocity)
    cost =  10 * abs(pi-final_angle) + 1000 * abs(final_velocity)
    return cost
end

# Monte Carlo simulation
function monte_carlo_optimization(num_simulations, initial_angle, max_time)
    best_cost = Inf
    best_solution = nothing
    
    degree = 5  # Polynomial degree for control function
    initial_coeffs = zeros(degree + 1)

    # Simulation loop
    for i in 1:num_simulations
        # Introduce random uncertainty in parameters
        mass = 1.0 + 0.1 * randn()   # Random mass perturbation ±10%
        length = 1.0 + 0.05 * randn() # Random length perturbation ±5%
        gravity = 9.81 + 0.1 * randn() # Random gravity perturbation ±1%
        
        # Initial conditions
        initial_conditions = [initial_angle, 0.0]  # Initial angle and velocity
       # parameters = [mass, length, gravity]

        parameters= coeffs, T 
        # Time span for the simulation
        tspan = (0.0, max_time)
        
        # Solve the system of equations
        prob = ODEProblem(pendulum_dynamics!, initial_conditions, tspan, parameters)
        solution = solve(prob, Tsit5(), saveat=0.01)
        
        # Evaluate the cost
        current_cost = cost_function(solution)
        
        # If the current cost is better, update the best solution
        if current_cost < best_cost
            best_cost = current_cost
            best_solution = solution
        end
    end
    
    return best_solution, best_cost
end

# Running the optimization
num_simulations = 1000
initial_angle = 0.0 # Start at 0.1 rad (5.7 degrees)
max_time = 10.0     # Max time for the simulation (seconds)

best_solution, best_cost = monte_carlo_optimization(num_simulations, initial_angle, max_time)

# Plot the best solution
if best_solution != nothing
    plot(best_solution.t, best_solution.u[1, :], label="Angle (θ)", xlabel="Time (s)", ylabel="Angle (rad)")
    plot!(best_solution.t, best_solution.u[2, :], label="Angular Velocity (θ_dot)", xlabel="Time (s)", ylabel="Velocity (rad/s)")
    title!("Optimized Inverted Pendulum Trajectory")
end
plot( best_solution, label="Angle (θ)", xlabel="Time (s)", ylabel="Angle (rad)")
println("Best Cost: ", best_cost)
