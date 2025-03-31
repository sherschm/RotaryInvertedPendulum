#WIP!
A_cl = sys_discrete.A - sys_discrete.B*K_opt  # Closed-loop A matrix

# Create the state-space system
sys_cl = ss(A_cl, sys_discrete.B, sys_discrete.C, 0)

sys_ol = sys_discrete

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