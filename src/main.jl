using RotaryInvertedPendulum
using LinearAlgebra
using DifferentialEquations
using Plots
using ControlSystems

#create model
include("parameters.jl")
(M,N,M_f,N_f), sys_cont, Total_energy, T_f,V_f = generate_dynamics(dyn_params,x_equil)

#calculate the Controller gains
Ts=0.01 #time step
sys_discrete=c2d(sys_cont,Ts)

#Define LQR objective parameters
Q=Diagonal([0.1,0.01,0.01,0.01]) 
R=0.001

#calculate LQR gains
K_opt=lqr(sys_discrete, Q, R) #Optimal feedback gain matrix

function u_f(x,t)
    # Acceleration control law using LQR
    acc = (-K_opt * (x - x_equil))[1]  # Control both angles
    # acc = (-K_opt * [0; (x - x_equil)[2:4]])[1]  # Control only θ2 stabilization

    # Apply acceleration limits
    max_acc = 200  # rad/s
    acc_limited = clamp(acc, -max_acc, max_acc)
    return acc_limited
end

# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl(x, p, t)
    θ=collect(x[1:2])
    θd=collect(x[3:4])
    
    Damping=p[1]
    D=[0 0;0 Damping]
    Damping_force=D*θd
  
    M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ...),N_f(x...),Damping_force)
    
    return vec([θd;inv(M_a)*(B_a*u_f(x,t)-N_a)])
end

#simulate!
q0=[0.1;pi+0.2;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ1d(t_0);θ2d(t_0)].
tspan = (0.0, 10.0)
ps=Damping
prob = ODEProblem(dynamics_acc_ctrl, q0, tspan,ps)

#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(prob)#,reltol=1e-10,abstol=1e-10)
plot(tvec,q_sol)
ctrl_vec=zeros(length(tvec))
energy=zeros(length(tvec))
for i in 1:length(tvec)
    ctrl_vec[i]=u_f([q_sol[i,:];qd_sol[i,:]],tvec[i])
    energy[i] = Total_energy([q_sol[i,:];qd_sol[i,:]])
end

#plot the free response of the generalised coordinates
plot(tvec,q_sol,label=["theta1" "theta2"],xlabel="Time (s)",ylabel="Angle (rad)")
savefig("plots//response")

plot(tvec,ctrl_vec,xlabel="Time (s)",ylabel="Control (rad/s^2)")
savefig("plots//control_acceleration")

#animate!
#plot_FFT(q_sol[:,2],tvec)
plot_params=(l1,l2)
rot_pendulum_animator(q_sol,tvec,plot_params;name="LQR_stabilisation")

##Calculate the swing-up trajectory!
q0=[0.0;0.0;0;0] #initial conditions
t_f=7
Δt=0.02 #trajectory time-step
n_traj=Int(round(t_f/Δt)) #number of trajectory points

tmax=n_traj*Δt
tvec=Δt:Δt:tmax
cmd=[x_equil[1];x_equil[2]] #pendulum command position
q_spin_up, qd_spin_up, qdd_spin_up,qddd_spin_up, torq_spin_up=SpinUpTrajectory(cmd,n_traj,Δt,q0,(M_f,N_f),T_f);

#animate!
rot_pendulum_animator(q_spin_up,tvec,plot_params;name="swing_up")
