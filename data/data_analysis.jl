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

include("../src/parameters.jl")

##symbolic model derivation - using Julia's Symbolics.jl toolbox
Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)  ud(t)#instantiate symbolic coordinates and their velocities (and accelerations), 
x=[θ1;θ2;θ1d;θ2d]

Dt = Differential(t) #time derivative operator

R20=Ry(θ2)*Rz(θ1) 

function correct_encoder_value(raw_value::Int, max_value::Int=65535)
    half_max = div(max_value + 1, 2) # Half of the max range
    return raw_value < half_max ? raw_value : raw_value - (max_value + 1)
end

function correct_encoder_vector(raw_values::Vector{Int}, max_value::Int=65535)
    return [correct_encoder_value(v, max_value) for v in raw_values]
end

#
data=CSV.read("C:\\Downloads\\500hz_SwingUp_Response.csv",DataFrame)
data=data[2:end,:]
#data=CSV.read("data/control_data_logs_example.csv",DataFrame)

#select nice data range
valid_i=6:size(data,1)-20
#valid_i=size(data,1)-20:size(data,1)
valid_i=1:size(data,1)

#t_delay=data.MessageRxTime[valid_i]-data.MessageTxTime[valid_i]

tvec=(data."#Time (uint32)[1]"[valid_i].-data."#Time (uint32)[1]"[valid_i[1]])/10^6

average_rate=length(tvec)/maximum(tvec) #Hz

θ1_raw=data."MotorPosition (int32)[1]"[valid_i]

θ2_raw=data."EncoderPosition (uint32)[1]"[valid_i]
encoder_points_per_rad=pi/1200

#convert to radians
θ2_vec=correct_encoder_vector(θ2_raw,65535)*encoder_points_per_rad
θ1_vec=θ1_raw*encoder_points_per_rad

plot(tvec,θ2_vec)

plot_rot_pendulum([θ1_vec[end];θ2_vec[end]],(l1,l2))

Δt=0.003
tvec_anim=minimum(tvec):0.05:maximum(tvec)

θ_data= general_lin_interp([θ1_vec θ2_vec],tvec,tvec_anim)

rot_pendulum_animator(θ_data,tvec_anim,(l1,l2);name="real_animation")

acc_ctrl_idxs=3000:3600

acc_cmd=data."RtAcc (float32)[1]"[valid_i]*encoder_points_per_rad
#plot(tvec[acc_ctrl_idxs],data.motor_Acceleration[acc_ctrl_idxs])
#plot(tvec[acc_ctrl_idxs],θ1_vec[acc_ctrl_idxs])

#numerical differentiate position data
window_size=9
poly_order=5

θ_ddot_diff=similar(θ_data)
for i in 1:2
    θ_ddot_diff[:,i]=savitzky_golay(θ_data[:,i], window_size, poly_order, deriv=2, rate=1/Δt).y
end

plot(tvec_anim,θ_ddot_diff[:,1],label="Real life",ylabel="Acceleration (encoder steps/s^2)",xlabel="time (s)")
acc_cmd_OG=readdlm("data/swing_up/swingup_acc_cmd.csv")
plot!(tvec,acc_cmd/encoder_points_per_rad,label="command_MARTe2")

plot!(acc_cmd_OG[:,1],acc_cmd_OG[:,2],label="command_ORIGINAL")
savefig("plots/swing_up_acc")


plot(tvec,acc_cmd/encoder_points_per_rad,label="command_MARTe2")
plot!(acc_cmd_OG[:,1],acc_cmd_OG[:,2],label="command_ORIGINAL")
savefig("plots/swing_up_acc_commands")

plot(tvec,θ2_vec,label="Real life",ylabel="Acceleration (encoder steps/s^2)",xlabel="time (s)")
pos_cmd=readdlm("data/swing_up/swingup_acc_cmd.csv")
#plot!(acc_cmd[:,1],acc_cmd[:,2],label="Simulation")
savefig("plots/swing_up_acc")


plot(tvec,acc_cmd,xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")


plot(p1,p2,layout=(2,1))

θ1_vec[acc_ctrl_idxs]


