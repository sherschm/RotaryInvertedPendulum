include("modelling_funcs.jl")
println("Creating model...")
include("parameters.jl")

##symbolic model derivation - using Julia's Symbolics.jl toolbox
Symbolics.@variables t θ1(t) θ2(t) θ1d(t) θ2d(t) θ1dd(t) θ2dd(t) u(t)  ud(t)#instantiate symbolic coordinates and their velocities (and accelerations), 
x=[θ1;θ2;θ1d;θ2d]

Dt = Differential(t) #time derivative operator

# rotation matrix of the body with respect to the inertial frame: a yaw-pitch rotation
R20=Ry(θ2)*Rz(θ1) 

#calculate the time derivative of R20 (will be useful in calculating the velocity of the centre of mass).
R20_dot=simplify.(expand_derivatives.(Dt.(R20,)))
R20_dot=substitute(R20_dot, Dict(Dt(θ1,) => θ1d, Dt(θ2,) => θ2d))

#calculate angular velocity of pendulum, with respect to inertial frame
Ω=R20_dot*R20'
ω=[Ω[3,2];Ω[1,3];Ω[2,1]]

#Rotational kinetic energy of pendulum L-bar
T_rot=0.5*ω'*Ip*ω

#Find the velocity of centre of mass of pendulum L-bar
v_com=R20_dot'*rc

#Linear Kinetic energy of pendulum L-bar centre of mass
T_lin=0.5*m*v_com'*v_com

#Total kinetic energy of pendulum L-bar
T=T_rot+T_lin

#Gravitational potential energy of pendulum L-bar
V=m*g*(R20'*rc)[3]

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

rot_pend_dynamics=rot_pend_dynamics_sym(ctrl_input_type,M,N,Damping)
x_equil= [0;pi;0;0]
ndof=2

#Calculate Linearised dynamics
A_sym = Symbolics.jacobian(rot_pend_dynamics, x)
subs_dict = Dict(vcat(x .=> x_equil, u => 0))
A = Float64.(Symbolics.substitute(A_sym, subs_dict)[1:2*ndof, 1:2*ndof])

B_sym = Symbolics.jacobian(rot_pend_dynamics, [u])
B = Float64.(Symbolics.substitute(B_sym, subs_dict))

C= I(4)#  #assume full observability [I(2) zeros(2,2)]
D=0

sys_cont=ss(A,B,C,D)