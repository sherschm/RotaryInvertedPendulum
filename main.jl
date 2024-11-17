
println("Importing packages...")
using Symbolics
using LinearAlgebra
using DifferentialEquations
using Plots
using Interpolations
using JuMP, Ipopt

include("plotting_funcs.jl")
include("modelling_funcs.jl")
include("TrajOptim.jl")

println("Creating model...")
##define system parameters
r=0.008 # L-bar radius (in m)
m1=0.008 # L-bar horizontal section mass (in kg)
m2=0.013 # L-bar vertical section length (in kg)
l1=0.08 # L-bar horizontal section length (in m)
l2=0.13 # L-bar vertical section length (in m)
m=m1+m2 #total L-bar mass
g=9.81 #gravity
rc=[0;0.06;-0.05] #position of L-bar centre of mass, measured in body frame (see diagrams)

#Inertia tensor of L-bar horizontal section measured in body frame (see diagrams)
I1=[(1/4)*m1*r^2+(1/3)*m1*l1^2 0 0;0 (1/2)*m1*r^2 0;0 0 (1/4)*m1*r^2+(1/3)m1*l1^2] 
#Inertia tensor of L-bar vertical section measured in body frame (see diagrams), using the parallel axis theorem
I2=[(1/4)*m2*r^2+(1/12)*m2*l2^2 0 0;0 (1/4)*m2*r^2+(1/12)*m2*l2^2 0;0 0 0.5*m2*r^2]+m2*([0;l1;-l2/2]'*[0;l1;-l2/2]*I(3)-[0;l1;-l2/2]*[0;l1;-l2/2]')
#Combined inertia tensor of L-bar, measured in body frame (see diagrams)
Ip=I1+I2

#Relevant rotation matrices
Rz(θ)=[cos(θ) -sin(θ) 0;
sin(θ) cos(θ) 0;
0 0 1] 

Ry(θ)=[cos(θ) 0 sin(θ);
0 1 0;
-sin(θ) 0 cos(θ)] 

##symbolic model derivation - using Julia's Symbolics.jl toolbox
Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) τ #instantiate symbolic coordinates and their velocities (and accelerations), 
Dt = Differential(t) #time derivative operator

# rotation matrix of the body with respect to the inertial frame: a yaw-pitch rotation
R20=Ry(θ2)*Rz(θ1) 

#calculate the time derivative of R20 (will be useful in calculating the velocity of the centre of mass).
R20_dot=simplify.(expand_derivatives.(Dt.(R20,)))
subs1=Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d)
R20_dot=substitute(R20_dot, (subs1))

#angular velocity of body, with respect to inertial frame
ω= Rz(θ1)'*[0;θ2d;0]+[0;0;θ1d]

#Rotational kinetic energy of pendulum L-bar
T_rot=0.5*ω'*Ip*ω

#velocity of centre of mass of pendulum L-bar
v_com=R20_dot'*rc

#Linear Kinetic energy of pendulum L-bar centre of mass
T_lin=0.5*m*v_com'*v_com

#Total kinetic energy of pendulum L-bar
T=T_rot+T_lin

#Gravitational potential energy of pendulum L-bar
V=m*g*(R20'*rc)[3]

#create functions for calculating system energy at a given time
T_f=eval(build_function(T, [θ1, θ2, θ1d, θ2d]...))
V_f=eval(build_function(V, [θ1, θ2, θ1d, θ2d]...))
Total_energy(x)=T_f(x...)+V_f(x...)

#generate Lagrangian model matrices / vectors
M,N=Lagrangian_dynamics(T,V)

##Convert symbolic terms M and N into Julia functions
M_f=eval(build_function(M, [θ1, θ2]...)[1])
N_f=eval(build_function(N, [θ1, θ2, θ1d, θ2d]...)[1])

#Define first order ODE form of dynamic model 
function rot_pend_dynamics(x,p,t)
    ## put dynamics in this form: xdot=f(x,u)
    τ=0 #e.g: no control
  
    Fric1=0# friction at the motor joint
    Fric2=0 # friction at the pendulum swing joint
  
    #Fric1=-0.001*x[3] # friction at the motor joint
    #Fric2=-0.00005*x[4] # friction at the pendulum swing joint
  
    return [x[3:4];inv(M_f(x[1:2]...))*([τ+Fric1;Fric2]-N_f(x...))]
end

#Define simulation problem parameters
q0=[0.0;0.1;0;0] #initial conditions
tspan = (0.0, 10.0)
prob = ODEProblem(rot_pend_dynamics, q0, tspan)

#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(prob)

#Plotting total energy can be a useful for model verification.
    #If no damping or actuation, energy should be constant,
    #   (or on the order of the ODE solver tolerance)
plot_energy(tvec,q_sol,qd_sol) 

#plot the response of the generalised coordinates
plot(tvec,q_sol,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Angle (rad)")
savefig("plots//response")

##Calculate the spin-up trajectory!
Δt=0.05 #trajectory time-step
n_traj=120 #number of trajectory points
tmax=n_traj*Δt
cmd=[0;pi] #pendulum command position
q_spin_up, qd_spin_up, qd_spin_up, torq_spin_up=SpinUpTrajectory(cmd,n_traj,Δt);





#=
## models of stepper motor system (e.g direct velocity control):
Symbolics.@variables u(t) ud(t)

M_stepper, N_stepper=stepper_dynamics(M,N)

M_stepper_f=eval(build_function(M_stepper, [θ1, θ2]))
N_stepper_f=eval(build_function(N_stepper, [θ1, θ2, u, θ2d, ud]))

#Define first order ODE form of dynamic model 
function rot_pend_dynamics_stepper(x,p,t)
    ## put dynamics in this form: xdot=f(x,u)

    Fric2=-0.00005*x[3] # friction at the pendulum swing joint
  
    return [input(t);x[4];inputd(t);(Fric2-N_stepper_f([x[1];x[2];input(t);x[4];inputd(t)]))/M_stepper_f([x[1];x[2]])]
end

input(t)=sin(10*t)
inputd(t)=10*cos(10*t)

q0=[0;0.1;input(0);0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, 10.0)
stepper_prob = ODEProblem(rot_pend_dynamics_stepper, q0, tspan)

#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(stepper_prob)

=#
## TO DO: ADD CONTROL SIMULATION - PID and then LQR, perhaps Energy shaping for 'spin up'?
## TO DO: find way to simulate forces from the stepper motor (since we can't control the torque directly)...
#  this may require Lagrange multipliers to constrain the velocity to a given value... not sure yet...

#=
#Define linearised system (linearised about unstable equilibrium, θ=pi )
A=
B=
C=
D= =#



