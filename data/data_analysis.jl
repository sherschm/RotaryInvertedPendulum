#script to read the output data from MARTe2 & animate the data

using DataFrames
using CSV
using Plots
using SavitzkyGolay
using LinearAlgebra
using Symbolics
using Interpolations

include("../src/parameters.jl")

#Relevant rotation matrices
Rz(θ)=[cos(θ) -sin(θ) 0;
sin(θ) cos(θ) 0;
0 0 1] 

Ry(θ)=[cos(θ) 0 sin(θ);
0 1 0;
-sin(θ) 0 cos(θ)] 

##symbolic model derivation - using Julia's Symbolics.jl toolbox
Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)  ud(t)#instantiate symbolic coordinates and their velocities (and accelerations), 
x=[θ1;θ2;θ1d;θ2d]

Dt = Differential(t) #time derivative operator

R20=Ry(θ2)*Rz(θ1) 

data=CSV.read("data/control_data_logs_example.csv",DataFrame)

#select nice data range
valid_i=6:size(data,1)-20
#valid_i=size(data,1)-20:size(data,1)
valid_i=1:size(data,1)

t_delay=data.MessageRxTime[valid_i]-data.MessageTxTime[valid_i]

tvec=(data.MessageRxTime[valid_i].-data.MessageRxTime[valid_i[1]])/10^9

average_rate=length(tvec)/maximum(tvec) #Hz

θ1_raw=data.rotor_position_steps[valid_i]

θ2_raw=data.encoder_position[valid_i]
encoder_points_per_rad=pi/1200

#convert to radians
θ2_vec=θ2_raw*encoder_points_per_rad
θ1_vec=θ1_raw*encoder_points_per_rad

plot(tvec,θ2_vec)

include("../src/plotting_funcs.jl")
plot_rot_pendulum([θ1_vec[end];θ2_vec[end]])

Δt=0.003
tvec_anim=minimum(tvec):0.05:maximum(tvec)

θ_data= general_lin_interp([θ1_vec θ2_vec],tvec,tvec_anim)

#rot_pendulum_animator(θ_data,tvec_anim;name="real_animation")

acc_ctrl_idxs=3000:3600
acc_cmd=data.motor_Acceleration[valid_i]*encoder_points_per_rad

plot(tvec[acc_ctrl_idxs],data.motor_Acceleration[acc_ctrl_idxs])
plot(tvec[acc_ctrl_idxs],θ1_vec[acc_ctrl_idxs])

#numerical differentiate position data
window_size=7
poly_order=5

θ_ddot_diff=similar(θ_data)
for i in 1:2
    θ_ddot_diff[:,i]=savitzky_golay(θ_data[:,i], window_size, poly_order, deriv=2, rate=1/Δt).y
end

p1=plot(tvec,acc_cmd)
p2=plot(tvec_anim,θ_ddot_diff[:,1])

plot(tvec[2300:2800],acc_cmd[2300:2800],xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")

plot(tvec,acc_cmd,xlabel="Time (s)",ylabel="Acceleration (rad/s^2)")


plot(p1,p2,layout=(2,1))

θ1_vec[acc_ctrl_idxs]
