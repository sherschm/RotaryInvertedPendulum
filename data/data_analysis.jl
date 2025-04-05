#script to read the output data from MARTe2 & animate the data

using RotaryInvertedPendulum
using DataFrames
using CSV
using Plots
using SavitzkyGolay
using LinearAlgebra
using Symbolics
using Interpolations
using DelimitedFiles
using DifferentialEquations

include("../src/parameters.jl")

function correct_encoder_value(raw_value::Int, max_value::Int=65535)
    half_max = div(max_value + 1, 2) # Half of the max range
    return raw_value < half_max ? raw_value : raw_value - (max_value + 1)
end

function correct_encoder_vector(raw_values::Vector{Int}, max_value::Int=65535)
    return [correct_encoder_value(v, max_value) for v in raw_values]
end

#
#data=CSV.read("C:\\Downloads\\test_result_200hz.csv",DataFrame)
data=CSV.read("C:\\Downloads\\test_result_reduced_jerk.csv",DataFrame)
data=data[2:end,:]
#data=CSV.read("data/control_data_logs_example.csv",DataFrame)

#select nice data range
valid_i=6:size(data,1)-20
#valid_i=size(data,1)-20:size(data,1)
valid_i=1:size(data,1)

#=
tvec=(data."#Time (uint32)[1]"[valid_i].-data."#Time (uint32)[1]"[valid_i[1]])/10^6

average_rate=length(tvec)/maximum(tvec) #Hz

θ1_raw=data."MotorPosition (int32)[1]"[valid_i]

θ2_raw=data."EncoderPosition (uint32)[1]"[valid_i]
encoder_points_per_rad=pi/1200

#convert to radians
#θ2_vec=correct_encoder_vector(θ2_raw,65535)*encoder_points_per_rad
θ2_vec=θ2_raw*encoder_points_per_rad.-θ2_raw[1]*encoder_points_per_rad
θ1_vec=θ1_raw*encoder_points_per_rad =#

θ1_vec=q_spin_up[:,1]
θ2_vec=q_spin_up[:,2]

plot(tvec,θ2_vec)

plot_rot_pendulum([θ1_vec[end];θ2_vec[end]],(l1,l2))

Δt=0.001

tvec_anim=minimum(tvec):Δt:maximum(tvec)

θ_data= general_lin_interp([θ1_vec θ2_vec],tvec,tvec_anim)

rot_pendulum_animator(θ_data,tvec_anim,(l1,l2);name="real_animation")

acc_cmd=data."RtAcc (float32)[1]"[valid_i]*encoder_points_per_rad

#numerical differentiate position data
window_size=11
poly_order=3

θ_ddot_diff=similar(θ_data)
for i in 1:2
    θ_ddot_diff[:,i]=savitzky_golay(θ_data[:,i], window_size, poly_order, deriv=2, rate=1/Δt).y
end

plot(tvec_anim,θ_ddot_diff[:,1],label="Real life",ylabel="Acceleration (encoder steps/s^2)",xlabel="time (s)")
acc_cmd_OG=readdlm("data/swing_up/swingup_acc_cmd.csv")
plot!(tvec,acc_cmd,label="command_MARTe2")

#plot!(acc_cmd_OG[:,1],acc_cmd_OG[:,2],label="command_ORIGINAL")
#savefig("plots/swing_up_acc")


plot(tvec,acc_cmd/encoder_points_per_rad,label="command_MARTe2")
#plot!(acc_cmd_OG[:,1],acc_cmd_OG[:,2],label="command_ORIGINAL")
#savefig("plots/swing_up_acc_commands")

plot(tvec,θ2_vec,label="Real life",ylabel="Acceleration (encoder steps/s^2)",xlabel="time (s)")
pos_cmd=readdlm("data/swing_up/swingup_acc_cmd.csv")
#plot!(acc_cmd[:,1],acc_cmd[:,2],label="Simulation")
savefig("plots/swing_up_acc")


plot(tvec,acc_cmd,xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")

#create model
(M,N,M_f,N_f), sys_cont, Total_energy, T_f,V_f = generate_dynamics(dyn_params,x_equil)


function interpolate_time_series(data::Matrix{Float64},time)
    N, k = size(data)
    #time = 1:N  # Assume uniform time steps from 1 to N
    interpolation_functions = [linear_interpolation(time, data[:, j]) for j in 1:k]

    # Return a function that interpolates for any time t
    return t -> [f(t) for f in interpolation_functions]
end

pos_cmd_f=interpolate_time_series(θ_data,tvec_anim)

#acc_cmd_f=interpolate_time_series(Matrix(θ_ddot_diff),tvec_anim)

acc_cmd_f=interpolate_time_series(θ_ddot_diff,tvec_anim)

# ODE function for DifferentialEquations.jl
function dynamics_acc_ctrl(x, p, t)
    θ=collect(x[1:2])
    θd=collect(x[3:4])
    
   acc_P_ctrl=100*(pos_cmd_f(t)[1]-θ[1])

    Damping=p[1]
    D=[0 0;0 Damping]
    Damping_force=D*θd
  
    M_a, N_a, B_a = dynamics_acc_ctrl_terms(M_f(θ...),N_f(x...),Damping_force)
    
    return vec([θd;inv(M_a)*(B_a*acc_P_ctrl-N_a)])
    #return vec([θd;inv(M_a)*(B_a*acc_cmd_f(t)[1]-N_a)])
end

#simulate!
q0=[0.0;0;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
#q0=[0.0;pi/10;0;0] #initial conditions - these are: [θ1(t_0);θ2(t_0);θ2d(t_0)].
tspan = (minimum(tvec), maximum(tvec))
ps=Damping
prob = ODEProblem(dynamics_acc_ctrl, q0, tspan,ps)

#Simulate & animate!
tvec,q_sol,qd_sol=pend_sim(prob)#,reltol=1e-10,abstol=1e-10)
plot(tvec,q_sol)
ctrl_vec=zeros(length(tvec))
energy=zeros(length(tvec))

plot_params=(l1,l2)
rot_pendulum_animator(q_sol,tvec,plot_params;name="sim_test")