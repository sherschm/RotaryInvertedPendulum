
#Relevant rotation matrices
Rz(θ)=[cos(θ) -sin(θ) 0;
sin(θ) cos(θ) 0;
0 0 1] 

Ry(θ)=[cos(θ) 0 sin(θ);
0 1 0;
-sin(θ) 0 cos(θ)] 

calc_R20(θ1,θ2)=Ry(θ2)*Rz(θ1) 

function Lagrangian_dynamics(T,V,vars)

  ##This calculates the terms of the equations of Motion & returns it in Symbolic form
  (t, θ1, θ2, θ1d, θ2d, θ1dd, θ2dd, u)=vars

  Symbolics.@variables t
  Dt = Differential(t) #time derivative operator

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
    subs=Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d,Dt(θ1d,) => θ1dd, Dt(θ2d,) => θ2dd)
    dt∂L∂qd[i]=substitute(dt∂L∂qd[i], (subs))
  end

  Eq=dt∂L∂qd-∂L∂q #Euler-lagrange equation without external control forces and damping.

  # create Mass matrix by separating out the acceleration terms to 
  M = Array{Num}(undef,2,2)
  M[1]=Symbolics.coeff(Eq[1], θ1dd)
  M[2]=Symbolics.coeff(Eq[1], θ2dd)
  M[3]=Symbolics.coeff(Eq[2], θ1dd)
  M[4]=Symbolics.coeff(Eq[2], θ2dd)

  # subtract mass matrix terms to find N vector (this is a vector describing the gravity and coriolis forces acting on/within the system)
  N=simplify.(expand.(simplify.(expand.(Eq-expand.(M*[θ1dd;θ2dd])))))

  return M, N
end


function pend_sim(prob)

  #Simulates and makes animation of the rotary pendulum system
  println("Generating response...")
  sol = solve(prob, reltol = 1e-10, abstol = 1e-10)#,Tsit5(), reltol=1e-8, abstol=1e-8)
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

function dynamics_acc_ctrl_terms(M,N,Damping_force)
  #Reduces the equations for the situation where velocity of θ1 is the control input.

  A=[1 0] #constraint

  M_acc=M
  N_acc=N+Damping_force-A'*inv(A*inv(M)*A')*A*inv(M)*(N+Damping_force)
  B_acc=A'*inv(A*inv(M)*A')

  return  M_acc, N_acc, B_acc
end

function rot_pend_dynamics_sym(ctrl_input_type,M,N,Damping)
  #This outputs the function that is integrated.
  #Slightly different depending on what the actual control input is

  Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)

  x=[θ1;θ2;θ1d;θ2d]
  D=[0 0;0 Damping]
  
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

    
    Damping_force=D*x[3:4]
    M_a,N_a,B_a=dynamics_acc_ctrl_terms(M,N,Damping_force)

    #Define first order ODE form of dynamic model 

    #Fric2=-0.00005*x[4] # friction at the pendulum swing joint
  
    #return [x[3];x[4];ctrl_law(t);(Fric2-N_stepper_f([x[1];x[2];x[3];x[4]]))/M_stepper_f([x[1];x[2]])]
    #return [x[3];x[4];ctrl_law(x,t);(Fric2-N_stepper/M_stepper)]
    return [x[3];x[4];inv(M_a)*(B_a*u-N_a)]
  elseif ctrl_input_type=="velocity"

    Symbolics.@variables u(t) ud(t)

    M_stepper, N_stepper=dynamics_vel_ctrl(M,N,u,ud)
    #M_stepper_f=eval(build_function(M_stepper, [θ1, θ2]))
    #N_stepper_f=eval(build_function(N_stepper, [θ1, θ2, u, θ2d, ud]))

    Fric2=-0.00005*x[3] # friction at the pendulum swing joint
  
    return [ctrl_law(x,t);x[4];ctrl_lawd(x,t);(Fric2-N_stepper_f([x[1];x[2];ctrl_law(x,t);x[4];ctrl_lawd(x,t)]))/M_stepper_f([x[1];x[2]])]

  end
end




function generate_dynamics(dyn_params,x_equil)
  rc,m,Ip,g,Damping=dyn_params

  ##symbolic model derivation - using Julia's Symbolics.jl toolbox
  Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)#instantiate symbolic coordinates and their velocities (and accelerations), 
  vars=(t, θ1, θ2, θ1d, θ2d, θ1dd, θ2dd, u)
  x=[θ1;θ2;θ1d;θ2d]

  Dt = Differential(t) #time derivative operator

  # rotation matrix of the body with respect to the inertial frame: a yaw-pitch rotation
  R20=calc_R20(θ1,θ2)

  #calculate the time derivative of R20 (will be useful in calculating the velocity of the centre of mass).
  R20_dot=simplify.(expand_derivatives.(Dt.(R20,)))
  R20_dot=substitute(R20_dot, Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d))

  #calculate angular velocity of pendulum, with respect to inertial frame
  Ω=R20_dot*R20'
  ω=[Ω[3,2];Ω[1,3];Ω[2,1]]

  #Rotational kinetic energy of pendulum L-bar
  T_rot=0.5*ω'*Ip*ω

  #Find the velocity of centre of mass of pendulum L-bar
  v_com=R20_dot'*rc

  #Linear Kinetic energy of pendulum L-bar centre of mass
  T_lin=0.5*m*v_com'*v_com

  #Total kinetic energy of pendulum L-bar
  T=T_rot+T_lin

  #Gravitational potential energy of pendulum L-bar
  V=m*g*(R20'*rc)[3]

  #create functions for calculating system energy at a given time
  T_f=eval(build_function(T, x...))
  V_f=eval(build_function(V, x...))
  Total_energy(x)=T_f(x...)+V_f(x...)

  #generate Lagrangian model matrices / vectors
  M,N=Lagrangian_dynamics(T,V,vars)

  ##Convert symbolic terms M and N into Julia functions
  M_f=eval(build_function(M, [θ1, θ2]...)[1])
  N_f=eval(build_function(N, x...)[1])

  ctrl_input_type="acceleration"#"torque"#

  rot_pend_dynamics=rot_pend_dynamics_sym(ctrl_input_type,M,N,Damping)
  ndof=2

  #Calculate Linearised dynamics
  A_sym = Symbolics.jacobian(rot_pend_dynamics, x)
  subs_dict = Dict(vcat(x .=> x_equil, u => 0))
  A = Float64.(Symbolics.substitute(A_sym, subs_dict)[1:2*ndof, 1:2*ndof])

  B_sym = Symbolics.jacobian(rot_pend_dynamics, [u])
  B = Float64.(Symbolics.substitute(B_sym, subs_dict))

  C= I(4)#  #assume full observability [I(2) zeros(2,2)]
  D=0

  linearised_sys_cont=ss(A,B,C,D)

  return (M,N,M_f,N_f), linearised_sys_cont, Total_energy, T_f,V_f 
end