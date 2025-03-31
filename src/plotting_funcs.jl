
using Statistics
using FFTW

L1_tip_pos(θ1,l1)=Rz(θ1)'*[0;l1;0]
L2_tip_pos(θ1,θ2,l1,l2)=calc_R20(θ1,θ2)'*[0;l1;-l2]

function plot_rot_pendulum(q,plot_params)

    l1,l2=plot_params
    coords=[zeros(3,1) L1_tip_pos(q[1],l1) L2_tip_pos(q[1],q[2],l1,l2)]

    #x_coords=[0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]]
    #y_coords=[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]]
    #z_coords=[0;L1_tip_pos(q[1])[3];L2_tip_pos(q[1],q[2])[3]]

    plot(coords[1,:],coords[2,:],coords[3,:],
    aspect_ratio=:equal,
    xlims=(-0.3,0.3),ylims=(-0.3,0.3),zlims=(-0.3,0.3),
    label=false,
    linewidth=8)
    #plot([0;L1_tip_pos(q[1])[1];L2_tip_pos(q[1],q[2])[1]],[0;L1_tip_pos(q[1])[2];L2_tip_pos(q[1],q[2])[2]],aspect_ratio=:equal,xlims=(-0.15,0.15),ylims=(-0.15,0.15))
end

function plot_rot_pendulum(q,i,Δt,plot_params)
    l1,l2=plot_params
    coords=[zeros(3,1) L1_tip_pos(q[1],l1) L2_tip_pos(q[1],q[2],l1,l2)]

    plot(coords[1,:],coords[2,:],coords[3,:],
    aspect_ratio=:equal,
    xlims=(-0.3,0.3),ylims=(-0.3,0.3),zlims=(-0.3,0.3),
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
  
function rot_pendulum_animator(x_sol,tvec,plot_params;name="rotary_pendulum_anim")
    #pendulum animation creation
    println("Creating animation...")
    anim_fps=20
    tvec_anim=1/anim_fps:1/anim_fps:maximum(tvec)

    x_anim=general_lin_interp(x_sol,tvec,tvec_anim)

    anim = @animate for i in 1:length(tvec_anim)
        plot()
        plot_rot_pendulum(x_anim[i,:],i,1/anim_fps,plot_params)
        plot!(label=string(i/anim_fps))
    end

    gif(anim,"anims//"*name*".gif",fps=anim_fps);
end


function plot_FFT(sig_in,tvec_in,file_name="spectrum")

  
    Δt_FFT=0.005
    Fs=1/Δt_FFT
    tvec_FFT=Δt_FFT:Δt_FFT:maximum(tvec_in)
    N=length(tvec_FFT)

    sig_interped=general_lin_interp(sig_in,tvec_in,tvec_FFT)

    sig=sig_interped.-mean(sig_interped)
    Y = fft(sig)
    #Y = fft(df.force_y)
    # Compute frequency vector
    frequencies = Fs * (0:(N ÷ 2 - 1)) / N

    # Compute the amplitude spectrum (single-sided)
    amplitude_spectrum = 2 * abs.(Y[1:N ÷ 2]) / N

    #find resonant peak
   #
   sorted_indices = sortperm(amplitude_spectrum, rev=true)
    println("fundamental_freq = "* string(frequencies[sorted_indices[1]]))

    fundamental_freq=frequencies[findfirst(item -> item == maximum(amplitude_spectrum),amplitude_spectrum)]

    # Plot the frequency spectrum
    plot(frequencies, amplitude_spectrum, xlabel="Frequency (Hz)", ylabel="Amplitude", title="Frequency Spectrum", legend=false,xlims=(0,10))
    #plot!([1.396;1.396],[0;maximum(amplitude_spectrum)],label="Observed Frequency")
    savefig("plots/"*file_name*".png")

    #writedlm("FFT_order-3.csv", [vec(frequencies) vec(amplitude_spectrum)], ',')
    return fundamental_freq 
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


#Plotting total energy can be a useful for model verification.
    #If no damping or actuation, energy should be constant,
    #   (or on the order of the ODE solver tolerance)
#plot_energy(tvec,q_sol,qd_sol) 