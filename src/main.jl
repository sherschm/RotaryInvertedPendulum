
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

#create model
include("modelling_script.jl")

Ts=0.001 #time step
sys_discrete=c2d(sys_cont,Ts)

Q=Diagonal([0.1,0.01,0.01,0.01]) 
R=0.01 #weighting matrix

K_opt=lqr(sys_discrete, Q, R) #Optimal feedback gain matrix
#K_opt=[-11.45 524.31 -16.75 89.38]
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
q0=[1;pi+0.2;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/10;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (0.0, 20.0)
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

plot_FFT(q_sol[:,2],tvec)
rot_pendulum_animator(q_sol,tvec;name="LQR_stabilisation")


##Calculate the spin-up trajectory!
q0=[0.0;0.0;0;0] #initial conditions
t_f=4
Δt=0.01 #trajectory time-step
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