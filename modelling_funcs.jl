
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
  sol = solve(prob, reltol = 1e-5, abstol = 1e-5)#,Tsit5(), reltol=1e-8, abstol=1e-8)
  tvec=sol.t

  #unpack solution
  x_sol=Array{Float64}(undef,length(tvec),2)  #joint position
  xdot_sol=Array{Float64}(undef,length(tvec),2) #joint velocity

  for i in 1:size(sol.u)[1]
      x_sol[i,:]=sol.u[i][1:2]
      xdot_sol[i,:]=sol.u[i][3:4]
  end

  #animate! 
  rot_pendulum_animator(x_sol,tvec) 
  return tvec, x_sol, xdot_sol
end


function stepper_dynamics(M,N)
  #Reduces the equations for the situation where velocity of θ1 is the control input.

  Eq=M[2,:]'*[θ1dd;θ2dd]+N[2]
  M_stepper=Symbolics.coeff(Eq, θ2dd)
  N_stepper=simplify.(expand.(Eq-expand.(M_stepper*θ2dd))) 

  subs=Dict(θ1d => u, θ1dd => ud)

  N_stepper=substitute(N_stepper,(subs))
  return M_stepper, N_stepper
end
