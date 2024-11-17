

L1_tip_pos(θ_1)=Rz(θ_1)'*[0;l1;0]
L2_tip_pos(θ_1,θ_2)=Symbolics.value.(substitute(R20,Dict([θ1 => θ_1, θ2 => θ_2])))'*[0;l1;-l2]


function plot_rot_pendulum(q)
    plot([0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]],[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]],[0;L1_tip_pos(q[1])[3];L2_tip_pos(q[1],q[2])[3]],
    aspect_ratio=:equal,
    xlims=(-0.15,0.15),ylims=(-0.15,0.15),zlims=(-0.15,0.15),
    label=false,
    linewidth=8)
    #plot([0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]],[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]],aspect_ratio=:equal,xlims=(-0.15,0.15),ylims=(-0.15,0.15))
end

function plot_rot_pendulum(q,i,Δt)
    plot([0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]],[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]],[0;L1_tip_pos(q[1])[3];L2_tip_pos(q[1],q[2])[3]],
    aspect_ratio=:equal,
    xlims=(-0.15,0.15),ylims=(-0.15,0.15),zlims=(-0.15,0.15),
    label=false,
    title="Time = "*string(floor(i*Δt))*" s",
    linewidth=8)
    #plot([0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]],[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]],aspect_ratio=:equal,xlims=(-0.15,0.15),ylims=(-0.15,0.15))
end

function general_lin_interp(dataset,tvec,tvec_new)
    ##Linear interpolation dataset from tvec to tvec_new
    #data input matrix should be N by X, where N is the time dimension. N can be any size and remains the same.

    interped_data=Array{Float64}(undef,length(tvec_new),size(dataset,2))
    for i in 1:size(dataset,2)
        interp_fn = LinearInterpolation(tvec, dataset[:,i],extrapolation_bc=Line()) 
        interped_data[:,i]=interp_fn(tvec_new)
    end
    return interped_data
end
  
function rot_pendulum_animator(x_sol,tvec;name="rotary_pendulum_anim")
    #pendulum animation creation
    println("Creating animation...")
    anim_fps=20
    tvec_anim=1/anim_fps:1/anim_fps:maximum(tvec)

    x_anim=general_lin_interp(x_sol,tvec,tvec_anim)

    anim = @animate for i in 1:length(tvec_anim)
        plot()
        plot_rot_pendulum(x_anim[i,:],i,1/anim_fps)
        plot!(label=string(i/anim_fps))
    end

    gif(anim,"anims//"*name*".gif",fps=anim_fps);
end

function plot_energy(tvec,q_sol,qd_sol)
    n=length(tvec)
    energy=zeros(n)
    for i in 1:n
        energy[i]=Total_energy([q_sol[i,:];qd_sol[i,:]])
    end
    plot(tvec,energy,xlabel="Time (s)",ylabel="Total system energy T+V (J)")
    savefig("plots//energy_plot")
end