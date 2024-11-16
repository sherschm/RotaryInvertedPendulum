using Symbolics
using LinearAlgebra
using Plots
using JuMP
using Ipopt
using DataFrames
using CSV

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

Δt=0.05
n_traj=120

cmd=[0.0;pi]

ndof=2

model = JuMP.Model(Ipopt.Optimizer)
#set_optimizer_attribute(model, "max_iter", 200)
#set_attribute(model, "linear_solver", "ma57")

#JuMP.@variable(model,q[i=1:ndof*n_traj],start=repeat(cmd,n_traj)[i])#,1:n_traj] )
JuMP.@variable(model,q[1:ndof*n_traj])#,1:n_traj] )
JuMP.@variable(model,qd[1:ndof*n_traj] )
JuMP.@variable(model,qdd[1:ndof*n_traj] )
JuMP.@variable(model,F_ee[1:n_traj] )

#=
function Static_obj_func(x_in::T...) where T

  q_in=collect(x_in[1:n_traj*ndof])
  qd_in=collect(x_in[n_traj*ndof.+(1:n_traj*ndof)])
  F_in=collect(x_in[2*n_traj*ndof.+(1:n_traj)])

  #KE_final=sum(qd_in[(j-1)*ndof+1:j*ndof]'*qd_in[(j-1)*ndof+1:j*ndof] for j in 1:n_traj)
  #T_err=sum(T_f(q_in[(j-1)*ndof+1:j*ndof]...,qd_in[(j-1)*ndof+1:j*ndof]...) for j in 1:n_traj)
  #V_err=-sum(V_f(q_in[(j-1)*ndof+1:j*ndof]...,qd_in[(j-1)*ndof+1:j*ndof]...) for j in 1:n_traj)
  tip_pos_err=sum((q_in[(j-1)*ndof+1:j*ndof][i]-cmd[i])^2 for i in 1:2, j in 1:n_traj)
  Objective = tip_pos_err #+KE_final  #+ tip_orientation_err+ torq_final#+ PE_final
  return T.(Objective)
end=#

function Static_obj_func(x_in::T...) where T

    q_in=collect(x_in[1:n_traj*ndof])
    qd_in=collect(x_in[n_traj*ndof.+(1:n_traj*ndof)])
    F_in=collect(x_in[2*n_traj*ndof.+(1:n_traj)])

    #KE_final=sum(qd_in[(j-1)*ndof+1:j*ndof]'*qd_in[(j-1)*ndof+1:j*ndof] for j in 1:n_traj)
    T_err=sum((Δt/2)*(T_f(q_in[(j-1)*ndof+1:j*ndof]...,qd_in[(j-1)*ndof+1:j*ndof]...)+T_f(q_in[(j)*ndof+1:(j+1)*ndof]...,qd_in[(j)*ndof+1:(j+1)*ndof]...)) for j in 1:n_traj-1)
    #V_err=-sum(V_f(q_in[(j-1)*ndof+1:j*ndof]...,qd_in[(j-1)*ndof+1:j*ndof]...) for j in 1:n_traj)
    #tip_pos_err=sum((q_in[(j-1)*ndof+1:j*ndof][i]-cmd[i])^2 for i in 1:2, j in 1:n_traj)
    Objective = T_err #+KE_final  #+ tip_orientation_err+ torq_final#+ PE_final
    return T.(Objective)
  end

function Dyn_constraint_vec(vars::T...) where T
    
  q_in=collect(vars[1:ndof])
  qd_in=collect(vars[ndof+1:2*ndof])
  qdd_in=collect(vars[2*ndof+1:3*ndof])

  F_in=collect(vars[3*ndof+1])

  dynamic_err = M_f(q_in...)*qdd_in +N_f(q_in...,qd_in...)-[F_in;0]#[3]

  return T.(dynamic_err)
end

Dyn_constraint_memo = memoize(Dyn_constraint_vec, ndof)

#JuMP.register(model, :objective, ndof+1, Static_obj_func; autodiff = true)
JuMP.@operator(model, objective, (2*ndof+1)*n_traj, Static_obj_func)
# JuMP.@operator(model, Dyn_constr, 2*ndof+6, Dyn_constraint)

#JuMP.@operator(model,  string_to_symbol_or_nothing("Dyn_constr1"), 2*ndof+6,Dyn_constraint_memo[1])
JuMP.@operator(model, Dyn_constr1, 3*ndof+1,Dyn_constraint_memo[1])
JuMP.@operator(model, Dyn_constr2, 3*ndof+1,Dyn_constraint_memo[2])

JuMP.@constraint(model,q[(n_traj-1)*ndof+1]==cmd[1])
JuMP.@constraint(model,q[(n_traj-1)*ndof+2]==cmd[2])
JuMP.@constraint(model,qd[(n_traj-1)*ndof+2]==0)
JuMP.@constraint(model,qd[(n_traj-1)*ndof+1]==0)

for j in 1:n_traj

  JuMP.@constraint(model,Dyn_constr1(q[(j-1)*ndof+1:j*ndof]...,qd[(j-1)*ndof+1:j*ndof]...,qdd[(j-1)*ndof+1:j*ndof]...,F_ee[j]...)==0)
  JuMP.@constraint(model,Dyn_constr2(q[(j-1)*ndof+1:j*ndof]...,qd[(j-1)*ndof+1:j*ndof]...,qdd[(j-1)*ndof+1:j*ndof]...,F_ee[j]...)==0)

  #end-effector velocity & acceleration constraints of stepper motor
  JuMP.@constraint(model,-0.4<=F_ee[j]<=0.4)
  #JuMP.@constraint(model,-pi<=q[(j-1)*ndof+1]<=pi)
 # JuMP.@constraint(model,-10<=qd[(j-1)*ndof+1]<=10)
  JuMP.@constraint(model,-10<=qdd[(j-1)*ndof+1]<=10)

  for i in 1:ndof
      if j==1
          JuMP.@constraint(model, q[i]==q0[i])
          JuMP.@constraint(model, qd[i,j]==q0[i+ndof])
      end
      if j>1
          JuMP.@constraint(model, 0.5*(qd[i+(j-1)*ndof]+qd[i+(j-2)*ndof])*Δt == q[i+(j-1)*ndof]-q[i+(j-2)*ndof])
          JuMP.@constraint(model, 0.5*(qdd[i+(j-1)*ndof]+qdd[i+(j-2)*ndof])*Δt == qd[i+(j-1)*ndof]-qd[i+(j-2)*ndof])
          # JuMP.@constraint(model, 0.5*(qd[i,j]+qd[i,j-1])*Δt == q[i,j]-q[i,j-1])
          #JuMP.@constraint(model, 0.5*(qdd[i,j]+qdd[i,j-1])*Δt == qd[i,j]-qd[i,j-1])
      end
      #JuMP.@constraint(model, -10<= q_eq_opt[i] <= 10)
  #JuMP.@constraint(model, q_eq_opt[i] == q_eq_real[i])
  end
end


JuMP.@objective(model, Min, objective(q...,qd...,F_ee...))
#JuMP.@objective(model, Min, 1)
#JuMP.@NLobjective(model, Min, Static_obj_func(Ehat,q_eq_opt...))
JuMP.optimize!(model)
#E_opt=JuMP.value.(Ehat)
if size(q,1)==1 || size(q,2)==1
    q_opt_val = reshape(JuMP.value.(q), (ndof, n_traj))'
    qd_opt_val = reshape(JuMP.value.(qd), (ndof, n_traj))'
    qdd_opt_val = reshape(JuMP.value.(qdd), (ndof, n_traj))'
    #F_opt_val=F_ee
else
    q_opt_val=JuMP.value.(q)
end

F_ee_opt=JuMP.value.(F_ee)

#JuMP.value.(J)
#return q_opt_val, JuMP.value.(F_ee)

rot_pendulum_animator(q_opt_val,Δt:Δt:n_traj*Δt)