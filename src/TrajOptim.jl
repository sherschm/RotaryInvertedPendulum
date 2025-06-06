using JuMP, Ipopt
using DelimitedFiles

function memoize(foo::Function, n_outputs::Int)
  last_x, last_f = nothing, nothing
  last_dx, last_dfdx = nothing, nothing
  function foo_i(i, x::T...) where {T<:Real}
      if T == Float64
          if x !== last_x
              last_x, last_f = x, foo(x...)
          end
          return last_f[i]::T
      else
          if x !== last_dx
              last_dx, last_dfdx = x, foo(x...)
          end
          return last_dfdx[i]::T
      end
  end
  return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end


function SpinUpTrajectory(cmd,n_traj,Δt,q0,dynamic_funcs,T_f)
    M_f, N_f=dynamic_funcs
    Damping=0.0002

    tvec=Δt:Δt:Δt*n_traj

    println("Building swing-up trajectory optimisation problem...")
    ndof=2
    model = JuMP.Model(Ipopt.Optimizer)
    #set_attribute(model, "linear_solver", "ma57")

    #Define optimisation variables. These are states at every trajectory time point
    JuMP.@variable(model,q[1:ndof*n_traj])
    JuMP.@variable(model,qd[1:ndof*n_traj] )
    JuMP.@variable(model,qdd[1:ndof*n_traj] )
    JuMP.@variable(model,qddd[1:n_traj] ) #motor jerk! limiting this smoothes the acceleration command

    JuMP.@variable(model,torq[1:n_traj] )

    function objective_func(x_in::T...) where T
        #Objective function: to minimise the kinetic energy throughout the trajectory
        q_in=collect(x_in[1:n_traj*ndof])
        qd_in=collect(x_in[n_traj*ndof.+(1:n_traj*ndof)])
        u_in=collect(x_in[2*n_traj*ndof.+(1:n_traj)])

        T_err=sum((Δt/2)*(T_f(q_in[(j-1)*ndof+1:j*ndof]...,qd_in[(j-1)*ndof+1:j*ndof]...)+T_f(q_in[(j)*ndof+1:(j+1)*ndof]...,qd_in[(j)*ndof+1:(j+1)*ndof]...)) for j in 1:n_traj-1)

        effort=sum(u_in[j]^2 for j in 1:n_traj-1)

        Objective = T_err#+effort#
        return T.(Objective)
    end
    
    function Dyn_constraint_vec(vars::T...) where T
        # A function that ensures the physics are considered in the optimisation
        # This will be used as a constraint.

        q_in=collect(vars[1:ndof])
        qd_in=collect(vars[ndof+1:2*ndof])
        qdd_in=collect(vars[2*ndof+1:3*ndof])
        F_in=collect(vars[3*ndof+1])

        D=[0 0;0 Damping] #damping matrix

        #This expression must equal zeros in the optimisation solution.
        dynamic_err = M_f(q_in...)*qdd_in +N_f(q_in...,qd_in...)-[F_in;0]-D*qd_in
     return T.(dynamic_err)
    end

#=
    # ODE function for DifferentialEquations.jl
    function Dyn_constraint_vec(vars::T...) where T
        θ_in=collect(vars[1:ndof])
        θd_in=collect(vars[ndof+1:2*ndof])
        θdd_in=collect(vars[2*ndof+1:3*ndof])
        u_in=collect(vars[3*ndof+1])

        x=[θ_in;θd_in]

        Damping=p[1]
        D=[0 0;0 Damping]
        Damping_force=D*θd_in
    
        M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ_in...),N_f(x...),Damping_force)
        
        return T.(M_a*θdd_in-(B_a*u_in-N_a))
    end=#

    #memoize the dynamic constraint. This is necessary in JuMP because Dyn_constraint_vec outputs a vector
    Dyn_constraint_memo = memoize(Dyn_constraint_vec, ndof)

    #Create nonlinear JuMP operators used in problem
    JuMP.@operator(model, objective, (2*ndof+1)*n_traj, objective_func)
    JuMP.@operator(model, Dyn_constr1, 3*ndof+1,Dyn_constraint_memo[1])
    JuMP.@operator(model, Dyn_constr2, 3*ndof+1,Dyn_constraint_memo[2])

    #constrain desired end-state
    JuMP.@constraint(model,q[(n_traj-1)*ndof+1]==cmd[1])
    JuMP.@constraint(model,q[(n_traj-1)*ndof+2]==cmd[2])
    JuMP.@constraint(model,qd[(n_traj-1)*ndof+2]==0)
    JuMP.@constraint(model,qd[(n_traj-1)*ndof+1]==0)

    #constrain start-state
    for i in 1:ndof
        JuMP.@constraint(model, q[i]==q0[i])
        JuMP.@constraint(model, qd[i]==q0[i+ndof])
    end

    for j in 1:n_traj
        #enforce dynamics 
        JuMP.@constraint(model,Dyn_constr1(q[(j-1)*ndof+1:j*ndof]...,qd[(j-1)*ndof+1:j*ndof]...,qdd[(j-1)*ndof+1:j*ndof]...,torq[j]...)==0)
        JuMP.@constraint(model,Dyn_constr2(q[(j-1)*ndof+1:j*ndof]...,qd[(j-1)*ndof+1:j*ndof]...,qdd[(j-1)*ndof+1:j*ndof]...,torq[j]...)==0)

        # define motor constraints if you need
        #STEPPER
        JuMP.@constraint(model,-0.4<=torq[j]<=0.4) 
        JuMP.@constraint(model,-20<=qdd[(j-1)*ndof+1]<=20) #acceleration constraints
        JuMP.@constraint(model,-200<=qddd[j]<=200) #jerk constraints

        if j>1
            JuMP.@constraint(model, 0.5*(qddd[j]+qddd[j-1])*Δt == qdd[1+(j-1)*ndof]-qdd[1+(j-2)*ndof])
            for i in 1:ndof
                #Trapezoidal rule
                JuMP.@constraint(model, 0.5*(qd[i+(j-1)*ndof]+qd[i+(j-2)*ndof])*Δt == q[i+(j-1)*ndof]-q[i+(j-2)*ndof])
                JuMP.@constraint(model, 0.5*(qdd[i+(j-1)*ndof]+qdd[i+(j-2)*ndof])*Δt == qd[i+(j-1)*ndof]-qd[i+(j-2)*ndof])
            end
        end
    end

    JuMP.@objective(model, Min, objective(q...,qd...,torq...))
    #Solve!
    JuMP.optimize!(model)

    #unpack results
    q_opt_val = reshape(JuMP.value.(q), (ndof, n_traj))'
    qd_opt_val = reshape(JuMP.value.(qd), (ndof, n_traj))'
    qdd_opt_val = reshape(JuMP.value.(qdd), (ndof, n_traj))'
    qddd_opt_val = JuMP.value.(qddd)
    torq_opt=JuMP.value.(torq)

    #plot the response of the generalised coordinates
    plot(tvec,q_opt_val,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Angle (rad)")
    savefig("plots//swing_up_plots//swing_up_traj")

    plot(tvec,torq_opt,xlabel="Time (s)",ylabel="Expected Torque (Nm)")
    savefig("plots//swing_up_plots//swing_up_torque")

    plot(tvec,qd_opt_val,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Velocity (rad/s)")
    savefig("plots//swing_up_plots//swing_up_velocity")

    plot(tvec,qdd_opt_val,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")
    savefig("plots//swing_up_plots//swing_up_accel")

    plot(tvec,qddd_opt_val,label=false,xlabel="Time (s)",ylabel="Motor Jerk (rad/s^3)")
    savefig("plots//swing_up_plots//swing_up_jerk")

    encoder_steps_per_rad=pi/1200

    writedlm("data/swing_up/swingup_pos.csv", [tvec q_opt_val],",")
    writedlm("data/swing_up/swingup_vel.csv", [tvec qd_opt_val],",")
    writedlm("data/swing_up/swingup_acc_cmd.csv", [Float32.(tvec) Float32.(qdd_opt_val[:,1]/encoder_steps_per_rad)],",")


    #output trajectory position, velocity, acceleration and motor torque profiles
    return q_opt_val, qd_opt_val, qdd_opt_val,qddd_opt_val, torq_opt
end
