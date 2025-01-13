function iLQR( x0, U, Q, R, Q_f; max_iters=100, tol=1e-6)
    N = length(U)               # Number of time steps
    n = length(x0)              # State dimension
    m = length(U[1])            # Control dimension

    X = [x0]                    # Initialize state trajectory
    for u in U                  # Forward pass (simulate trajectory)
        push!(X, dynamics(X[end], u))
    end

    for iter in 1:max_iters
        # Backward pass
        K = Vector{Matrix}(undef, N)  # Feedback gains
        k = Vector{Vector}(undef, N) # Feedforward terms

        P = Q_f  # Terminal cost-to-go
        V_x = Q_f * X[end]
        V_xx = Q_f

        for k in N:-1:1
            xk, uk = X[k], U[k]

            A=A_f(xk)
            B=B_f(xk)

            Q_x = Q * xk
            Q_u = R * uk
            Q_xx = Q
            Q_ux = zeros(m, n)
            Q_uu = R

            K[k] = -inv(Q_uu) * Q_ux
            k[k] = -inv(Q_uu) * Q_u
            P = Q_xx + A' * P * A + K[k]' * Q_uu * K[k]
        end

        # Forward pass to update U and X
        new_U = []
        new_X = [x0]
        for k in 1:N
            uk = U[k] + k[k] + K[k] * (new_X[end] - X[k])
            push!(new_U, uk)
            push!(new_X, 
            
            new_X[end]+Δt_ctrl*ODE([new_X[end];0]) dynamics(, uk))
        end

        # Check convergence
        if norm(U - new_U) < tol
            break
        end
        U = new_U
        X = new_X
    end

    return X, U
end