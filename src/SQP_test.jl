
using LinearAlgebra, SparseArrays, OSQP

function compute_linearized_dynamics(f, g, ∂f, ∂g, xk, uk, Δt)
    # Evaluate f(x) and g(x) at the current state
    fx = f(xk...)
    gx = g(xk...)

    # Compute the Jacobians
    ∂fx = ∂f(xk...)         # Jacobian of f(x) w.r.t. x
    ∂gx = ∂g(xk...)         # Jacobian of g(x) w.r.t. x

    # Compute A_k and B_k
    Ak = ∂fx + ∂gx * uk     # Linearized state dynamics w.r.t x
    Bk = gx                 # Linearized control dynamics

    # Nominal dynamics (used for the affine term)
    nominal_dynamics = fx + gx * uk

    # Discretize the dynamics
    Ad = I + Δt * Ak        # Discrete-time state matrix
    Bd = Δt * Bk            # Discrete-time control matrix
    bd = Δt * nominal_dynamics  # Affine term for discretized dynamics

    return Ad, Bd, bd
end

function solve_mpc_osqp(f, g, ∂f, ∂g, x0, x_ref, u_bounds,vel_bounds, Q, R, Δt, N)

    nx = length(x0)       # Number of states
    nu = 1  # Number of control inputs

    # Total decision variables: x and u over the horizon
    n_decision_vars = N * (nx + nu)

    # Initialize QP matrices
    P = spzeros(n_decision_vars, n_decision_vars)  # Quadratic cost
    q = zeros(n_decision_vars)                    # Linear cost
    A_dyn = spzeros(N * nx, n_decision_vars)          # Dynamics constraints
    lb_dyn = zeros(N * nx)                             # Lower bounds
    ub_dyn = zeros(N * nx)                             # Upper bounds

    A_torq = spzeros(N * nu, n_decision_vars)          # Dynamics constraints

    # Initialize state and control variables
    xk = x0
    idx_x = 1
    idx_u = nx + 1

    for k in 1:N
        # Compute linearized dynamics at (xk, uk)
        Ad, Bd, bd = compute_linearized_dynamics(f, g, ∂f, ∂g, xk, 0.0, Δt)

        # Dynamics constraints for timestep k
        A_dyn[(k-1)*nx+1:k*nx, idx_x:idx_x+nx-1] = -(I + Δt * Ad)
        A_dyn[(k-1)*nx+1:k*nx, idx_u:idx_u+nu-1] = -Δt * Bd
        if k<N
            A_dyn[(k-1)*nx+1:k*nx, idx_x+nx:idx_x+2*nx-1] = I(nx)
        end
        lb_dyn[(k-1)*nx+1:k*nx] = -bd
        ub_dyn[(k-1)*nx+1:k*nx] = -bd

        # Cost function: tracking (x - x_ref)
        Q_block = Δt * Q
        R_block = Δt * R
        P[idx_x:idx_x+nx-1, idx_x:idx_x+nx-1] .= Q_block
        P[idx_u:idx_u+nu-1, idx_u:idx_u+nu-1] .= R_block

        #P[idx_x:idx_x+nx-1, idx_x:idx_x+nx-1] .+ Δt *[I(2) zeros(2,2);zeros(2,4)] 
        #q[idx_x:idx_x+nx-1] .= -2 * Δt * [I(2) zeros(2,2);zeros(2,4)] * x_ref  # Tracking term
        q[idx_x:idx_x+nx-1] .= -2 * Δt * Q * x_ref  # Tracking term

        A_torq[k,idx_u:idx_u+nu-1].=1.0;

        # Update indices for next timestep
        idx_x += nx + nu
        idx_u += nx + nu

    end

    A_IC=spzeros(nx, n_decision_vars)
    A_IC[1:nx,1:nx]=I(nx)
    lb_IC=x0
    ub_IC=x0

    lb_torq=vec(repeat([-u_bounds],N*nu,1))
    ub_torq=vec(repeat([u_bounds],N*nu,1))

    A=[A_dyn;A_IC;A_torq]
    lb=[lb_dyn;lb_IC;lb_torq]
    ub=[ub_dyn;ub_IC;ub_torq]


    # OSQP setup
    prob = OSQP.Model()
    OSQP.setup!(prob, P=P, q=q, A=A, l=lb, u=ub, polish=true)

    # Solve the QP
    result = OSQP.solve!(prob)

    # Extract solution
    z_opt = result.x
    
    idx_x = 1
    idx_u = nx + 1
    x_opt=Array{Float64}(undef,N,nx)
    u_opt=Array{Float64}(undef,N,nu)

    for k in 1:N
        x_opt[k,:] = z_opt[idx_x:idx_x+nx-1]
        u_opt[k,:] = z_opt[idx_u:idx_u+nu-1]
        idx_x += nx + nu
        idx_u += nx + nu
    end

    return x_opt, u_opt
end


# Horizon parameters
N =20    # Horizon length
Δt = 0.001            # Time step
ndof = 2             # Degrees of freedom

x0=[0.0;0.0;0.0;0.0]
x_ref=[0.0;pi;0.0;0.0]

# Bounds
τ_min = -0.4; τ_max = 0.4
q_min = -Inf; q_max = Inf
qd_min = -10.0; qd_max = 10.0

u_bounds=0.4
vel_bounds=2.0
R=0.0
#Q=[zeros(2,4);zeros(2,2) I(2)]
Q= I(4)

