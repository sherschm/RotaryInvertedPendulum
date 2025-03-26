
println("Importing packages...")

using RotaryInvertedPendulum
using Symbolics
using LinearAlgebra
using DifferentialEquations
using Plots
using Interpolations
using ControlSystems
using DelimitedFiles

include("plotting_funcs.jl")
include("modelling_funcs.jl")

println("Creating model...")
include("parameters.jl")

##symbolic model derivation - using Julia's Symbolics.jl toolbox
Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)  ud(t)#instantiate symbolic coordinates and their velocities (and accelerations), 
x=[θ1;θ2;θ1d;θ2d]

Dt = Differential(t) #time derivative operator

# rotation matrix of the body with respect to the inertial frame: a yaw-pitch rotation
R02=Rz(θ1)*Ry(θ2)

#calculate the time derivative of R20 (will be useful in calculating the velocity of the centre of mass).
R02_dot=simplify.(expand_derivatives.(Dt.(R02,)))
subs1=Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d)
R02_dot=substitute(R02_dot, (subs1))

#angular velocity of body, with respect to inertial frame
#ω= Rz(θ1)'*[0;θ2d;0]+[0;0;θ1d]
#ω=[(R20'*R20_dot)[3,2];(R20'*R20_dot)[1,3];(R20'*R20_dot)[2,1]]
ω=[(R02_dot*R02')[3,2];(R02_dot*R02')[1,3];(R02_dot*R02')[2,1]]

#Rotational kinetic energy of pendulum L-bar
T_rot=0.5*ω'*Ip*ω

#velocity of centre of mass of pendulum L-bar
v_com=R02_dot*rc

#v_com=cross(ω,rc)
v_com_f=eval(build_function(v_com, x...)[1])

#Linear Kinetic energy of pendulum L-bar centre of mass
T_lin=0.5*m*v_com'*v_com

#Total kinetic energy of pendulum L-bar
T=T_rot+T_lin

#Gravitational potential energy of pendulum L-bar
V=m*g*(R02*rc)[3]

#create functions for calculating system energy at a given time
T_f=eval(build_function(T, x...))
V_f=eval(build_function(V, x...))
Total_energy(x)=T_f(x...)+V_f(x...)

#generate Lagrangian model matrices / vectors
M,N=Lagrangian_dynamics(T,V)

##Convert symbolic terms M and N into Julia functions
M_f=eval(build_function(M, [θ1, θ2]...)[1])
N_f=eval(build_function(N, x...)[1])
dynamic_funcs=(M_f,N_f)

M_jac=Symbolics.jacobian(M,[θ1, θ2])
N_jac=Symbolics.jacobian(N,x)

ctrl_input_type="acceleration"#"torque"#

Damping = [0 0;0 0.0001]

rot_pend_dynamics=rot_pend_dynamics_sym(ctrl_input_type,M,N,Damping)
x_equil= [0;pi;0;0]
ndof=2

##SIMULATE LQR
#Calculate Linearised dynamics
A_f=eval(build_function(Symbolics.jacobian(rot_pend_dynamics,x),[x;u])[1])
B_f=eval(build_function(Symbolics.jacobian(rot_pend_dynamics,[u]),x)[1])

A=A_f([x_equil;0])[1:2*ndof,1:2*ndof]
B=B_f(x_equil)
C= [I(2) zeros(2,2)] #assume full observabilityI(4)
D=0
sys_cont=ss(A,B,C,D)


Ts=0.001
sys_discrete=c2d(sys_cont,Ts)

Q=Diagonal([0.1,0.01,0.01,0.01]) 
R=0.01 #weighting matrix
#K_opt=lqr(sys_cont, Q, R) #Optimal feedback gain matrix
K_opt=lqr(sys_discrete, Q, R) #Optimal feedback gain matrix

#K_opt=[ -0.000994769  -1010.51  -0.318096  -133.896]
#K_opt=[11.45 524.31 16.75 89.38]
#K_opt=[-0.15 500 -4 100]
#K_opt=[0 6 0 4.8]
#K_opt=[-0.5 300 -10 30]
A_cl = A - B*K_opt  # Closed-loop A matrix
C = I(4)#[1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0]  #     # Output matrix (adjust based on your system)
#D = [0]         # Direct transmission term

# Create the state-space system
sys_cl = ss(A_cl, B, C, 0)

sys_ol = ss(A, B, C, 0)

# Convert to transfer function (optional, useful for visualization)

sys_ol_tf = tf(sys_ol)
ol_poles=poles(sys_ol_tf)

sys_tf = tf(sys_cl)

cl_poles=poles(sys_tf)
w = exp10.(LinRange(-8, 7, 1000))  # Log-spaced frequency range
p_nyq1=nyquistplot(sys_tf[1,1],w,linewidth=7,title="theta1")#,xlims=(-2,2),ylims=(-2.5,2))
p_nyq2=nyquistplot(sys_tf[2,1],w,linewidth=7,title="theta2")
p_nyq3=nyquistplot(sys_tf[3,1],w,linewidth=7,title="theta1d")
p_nyq4=nyquistplot(sys_tf[4,1],w,linewidth=7,title="theta2d")

plot(p_nyq1,p_nyq2,p_nyq3,p_nyq4,layout=(2,2))
savefig("plots/Nyquist")

bodeplot(sys_tf[3,1])
plot(step(sys_tf[4,1]))

marginplot(sys_tf[1,1])

#K_opt=[0 300 0 30]

function u_f(x,t)
    ##acceleration control law

    #LQR:
    acc=(-K_opt*(x-x_equil))[1] #if we want to control both angles
   # acc=(-K_opt*([0;(x-x_equil)[2:4]]))[1] #if we only care about θ2 stabilisation

    #apply maximum acceleration limits
    max_acc=200 #rad/s

    if acc>=max_acc
        acc_limited=max_acc
    elseif acc<=-max_acc
        acc_limited=-max_acc
    else
        acc_limited=acc
    end

    #return acc
    return acc_limited
    #return 0
end

#simulate!
q0=[-0.0;pi+0.01;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/8;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, 10.0)
prob = ODEProblem(dynamics_acc_ctrl, q0, tspan)

#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(prob)#,reltol=1e-10,abstol=1e-10)
plot(tvec,q_sol)
ctrl_vec=zeros(length(tvec))
energy=zeros(length(tvec))
for i in 1:length(tvec)
    ctrl_vec[i]=u_f([q_sol[i,:];qd_sol[i,:]],tvec[i])
    energy[i] = Total_energy([q_sol[i,:];qd_sol[i,:]])
end

#plot the response of the generalised coordinates
plot(tvec,q_sol,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Angle (rad)")
savefig("plots//response")

plot(tvec,ctrl_vec,xlabel="Time (s)",ylabel="Control (rad/s^2)")
savefig("plots//control_acceleration")

#animate!
rot_pendulum_animator(q_sol,tvec;name="LQR_stabilisation")

##Calculate the spin-up trajectory!
q0=[0.0;0.0;0;0] #initial conditions
t_f=2
Δt=0.05 #trajectory time-step
n_traj=Int(round(t_f/Δt)) #number of trajectory points

tmax=n_traj*Δt
tvec=Δt:Δt:tmax
cmd=[x_equil[1];x_equil[2]] #pendulum command position
q_spin_up, qd_spin_up, qdd_spin_up,qddd_spin_up, torq_spin_up=SpinUpTrajectory(cmd,n_traj,Δt,q0,dynamic_funcs,T_f);

#animate!
rot_pendulum_animator(q_spin_up,tvec;name="swing_up")

#plot the response of the generalised coordinates
plot(tvec,q_spin_up,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Angle (rad)")
savefig("plots//swing_up_plots//swing_up_traj")

plot(tvec,torq_spin_up,xlabel="Time (s)",ylabel="Expected Torque (Nm)")
savefig("plots//swing_up_plots//swing_up_torque")

plot(tvec,qd_spin_up,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Velocity (rad/s)")
savefig("plots//swing_up_plots//swing_up_velocity")

plot(tvec,qdd_spin_up,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")
savefig("plots//swing_up_plots//swing_up_accel")

plot(tvec,qddd_spin_up,label=false,xlabel="Time (s)",ylabel="Motor Jerk (rad/s^3)")
savefig("plots//swing_up_plots//swing_up_jerk")

encoder_points_per_rad=pi/1200

writedlm("data/swing_up/swingup_pos_cmd.csv", [tvec q_spin_up[:,1]])
writedlm("data/swing_up/swingup_vel_cmd.csv", [tvec qd_spin_up[:,1]])
writedlm("data/swing_up/swingup_acc_cmd.csv", [Float32.(tvec) Float32.(qdd_spin_up[:,1]/encoder_points_per_rad)])