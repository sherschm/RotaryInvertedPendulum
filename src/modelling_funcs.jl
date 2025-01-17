
#Relevant rotation matrices
Rz(θ)=[cos(θ) -sin(θ) 0;
sin(θ) cos(θ) 0;
0 0 1] 

Ry(θ)=[cos(θ) 0 sin(θ);
0 1 0;
-sin(θ) 0 cos(θ)] 


function Lagrangian_dynamics(T,V)

  ##This calculates the terms of the equations of Motion & returns it in Symbolic form

  #System Lagrangian
  L=T-V
  #Derive Symbolic Lagrangian Equations of Motion 
  Dθ1 = Differential(θ1) #derivative operator
  Dθ2 = Differential(θ2) #derivative operator
  Dθ1d = Differential(θ1d) #derivative operator
  Dθ2d = Differential(θ2d) #derivative operator
  ∂L∂q= simplify.(expand.(expand_derivatives.([Dθ1(L);Dθ2(L)])))
  ∂L∂qd= simplify.(expand.(expand_derivatives.([Dθ1d(L);Dθ2d(L)])))
  dt∂L∂qd=simplify.(expand.(expand_derivatives.(Dt.(∂L∂qd))))

  for i in 1:length(dt∂L∂qd)
    subs2=Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d,Dt(θ1d,) => θ1dd, Dt(θ2d,) => θ2dd)
    dt∂L∂qd[i]=substitute(dt∂L∂qd[i], (subs2))
  end

  Eq=dt∂L∂qd-∂L∂q #Euler-lagrange equation without external control forces and damping.

  # create Mass matrix by separating out the acceleration terms to 
  M = Array{Num}(undef,2,2)
  M[1]=Symbolics.coeff(Eq[1], θ1dd)
  M[2]=Symbolics.coeff(Eq[1], θ2dd)
  M[3]=Symbolics.coeff(Eq[2], θ1dd)
  M[4]=Symbolics.coeff(Eq[2], θ2dd)

  # subtract mass matrix terms to find N vector (this is a vector describing the gravity and coriolis forces acting on/within the system)
  N=simplify.(expand.(Eq-expand.(M*[θ1dd;θ2dd]))) 
  return M, N
end


function pend_sim(prob)

  #Simulates and makes animation of the rotary pendulum system
  println("Generating response...")
  sol = solve(prob)#, reltol = 1e-5, abstol = 1e-5)#,Tsit5(), reltol=1e-8, abstol=1e-8)
  tvec=sol.t

  #unpack solution
  x_sol=Array{Float64}(undef,length(tvec),2)  #joint position
  xdot_sol=Array{Float64}(undef,length(tvec),2) #joint velocity

  for i in 1:size(sol.u)[1]
      x_sol[i,:]=sol.u[i][1:2]
      xdot_sol[i,:]=sol.u[i][3:4]
  end

  #animate! 
  return tvec, x_sol, xdot_sol
end


function dynamics_vel_ctrl(M,N,u,ud)
  #Reduces the equations for the situation where velocity of θ1 is the control input.

  Eq=M[2,:]'*[θ1dd;θ2dd]+N[2]
  M_stepper=Symbolics.coeff(Eq, θ2dd)
  N_stepper=simplify.(expand.(Eq-expand.(M_stepper*θ2dd))) 

  subs=Dict(θ1d => u, θ1dd => ud)

  N_stepper=substitute(N_stepper,(subs))
  return M_stepper, N_stepper
end

function dynamics_acc_ctrl_terms(M,N)
  #Reduces the equations for the situation where velocity of θ1 is the control input.

  A=[1 0] #constraint

  M_acc=M
  N_acc=N-A'*inv(A*inv(M)*A')*A*inv(M)*N
  B_acc=A'*inv(A*inv(M)*A')

  return  M_acc, N_acc, B_acc
end

function rot_pend_dynamics_sym(ctrl_input_type,M,N)
  #This outputs the function that is integrated.
  #Slightly different depending on what the actual control input is

  Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)

  x=[θ1;θ2;θ1d;θ2d]

  if ctrl_input_type=="torque"

    ## put dynamics in this form: xdot=f(x,u)
    τ=0 #e.g: no control

    Fric1=0# friction at the motor joint
    Fric2=0 # friction at the pendulum swing joint

    #Fric1=-0.001*x[3] # friction at the motor joint
    #Fric2=-0.00005*x[4] # friction at the pendulum swing joint

    return [x[3:4];inv(M)*([u+Fric1;Fric2]-N)]
   #return [x[3:4];inv(M(x[1:2]...))*([u+Fric1;Fric2]-N(x...))]

  elseif ctrl_input_type=="acceleration"
    
    M_a,N_a,B_a=dynamics_acc_ctrl_terms(M,N)

    #Define first order ODE form of dynamic model 

    Fric2=0.0#-0.00005*x[3] # friction at the pendulum swing joint
  
    #return [x[3];x[4];ctrl_law(t);(Fric2-N_stepper_f([x[1];x[2];x[3];x[4]]))/M_stepper_f([x[1];x[2]])]
    #return [x[3];x[4];ctrl_law(x,t);(Fric2-N_stepper/M_stepper)]
    return [x[3];x[4];inv(M_a)*(B_a*u-[0.0;Fric2]-N_a)]
  elseif ctrl_input_type=="velocity"

    Symbolics.@variables u(t) ud(t)

    M_stepper, N_stepper=dynamics_vel_ctrl(M,N,u,ud)
    #M_stepper_f=eval(build_function(M_stepper, [θ1, θ2]))
    #N_stepper_f=eval(build_function(N_stepper, [θ1, θ2, u, θ2d, ud]))

    Fric2=-0.00005*x[3] # friction at the pendulum swing joint
  
    return [ctrl_law(x,t);x[4];ctrl_lawd(x,t);(Fric2-N_stepper_f([x[1];x[2];ctrl_law(x,t);x[4];ctrl_lawd(x,t)]))/M_stepper_f([x[1];x[2]])]

  end
end


# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl(x, p, t)
  θ=collect(x[1:2])
  θd=collect(x[3:4])
  M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ...),N_f(x...))


  Fric2=0.0#-0.00005*x[3] # friction at the pendulum swing joint

  return vec([θd;inv(M_a)*(B_a*u_f(x,t)-[0.0;Fric2]-N_a)])
end