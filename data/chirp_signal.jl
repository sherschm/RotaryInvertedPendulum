
using ChirpSignal
using DelimitedFiles
using Plots
using CSV

#trajectory parameters
max_freq=15#2*pi*5 #omega_max
tmax=50#s
fs=1000 # signal sample rate
f0=0.1  # signal starting frequency
f1=12  # signal final frequency
Δt=1/fs
tvec=Δt:Δt:tmax
n=length(tvec)

#generate velocity trajectory
vel_cmd=10*chirp(tmax, fs, f0, f1; method="linear", phase=0.0)

#calculate position traj
pos=similar(vel_cmd)
for i in 2:n
    pos[i]=pos[i-1]+Δt*vel_cmd[i]
end

#calculate position traj
acc=similar(vel_cmd)
for i in 2:n
    acc[i]=(vel_cmd[i]-vel_cmd[i-1])/Δt
end

#plot trajectories
p_pos=plot(tvec,pos,ylabel="Position (rad)",label=false)
p_vel=plot(tvec,vel_cmd,ylabel="Velocity (rad/s)",label=false)
p_acc=plot(tvec,acc,xlabel="Time (s)",ylabel="Acceleration (rad/s^2)",label=false)

plot(p_pos,p_vel,p_acc,layout=(3,1))

#Output trajectory to csv
Motor_pos=[tvec (1200/pi)*pos]
Motor_vel=[tvec (1200/pi)*vel_cmd]
Motor_acc=[tvec (1200/pi)*acc]

writedlm("data/chirp_position_cmd.csv", Motor_pos)
writedlm("data/chirp_velocity_cmd.csv", Motor_vel)
writedlm("data/chirp_acceleration_cmd.csv", Motor_acc)

p_pos=plot(Motor_pos[:,1],Motor_pos[:,2],xlabel="Time (s)",ylabel="Position (steps)",linewidth=2,label=false)
p_vel=plot(Motor_vel[:,1],Motor_vel[:,2],xlabel="Time (s)",ylabel="Velocity (steps/s)",linewidth=2,label=false)
p_acc=plot(Motor_acc[:,1],Motor_acc[:,2],xlabel="Time (s)",ylabel="Acceleration (steps/s^2)",linewidth=2,label=false)

plot(p_pos,p_vel,p_acc,layout=(3,1),size=(1000,1000),dpi=300)
savefig("data/chirp_cmd")