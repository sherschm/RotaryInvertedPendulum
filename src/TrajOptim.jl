using JuMP, Ipopt

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

        #This expression must equal zeros in the optimisation solution.
        dynamic_err = M_f(q_in...)*qdd_in +N_f(q_in...,qd_in...)-[F_in;0]#+damping if used

     return T.(dynamic_err)
    end

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

        # define motor velocity & acceleration constraints if you need
        #STEPPER
        JuMP.@constraint(model,-0.4<=torq[j]<=0.4)
        JuMP.@constraint(model,-20<=qdd[(j-1)*ndof+1]<=20)

        JuMP.@constraint(model,-500<=qddd[j]<=500)
        #Dynamixel Motor XM540-W270-T/R
       # JuMP.@constraint(model,-10.6<=torq[j]<=10.6)
        #JuMP.@constraint(model,-10<=qdd[(j-1)*ndof+1]<=10)
        #JuMP.@constraint(model,-pi<=qd[(j-1)*ndof+1]<=pi)
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

    #tvec=Δt:Δt:n_traj*Δt

    #output trajectory position, velocity, acceleration and motor torque profiles
    return q_opt_val, qd_opt_val, qdd_opt_val,qddd_opt_val, torq_opt
end
