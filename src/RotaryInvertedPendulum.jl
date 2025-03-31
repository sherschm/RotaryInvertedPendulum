__precompile__()
module RotaryInvertedPendulum

    using Symbolics
    using LinearAlgebra
    using DifferentialEquations
    using Plots
    using Interpolations
    using ControlSystems
    using JuMP
    using Ipopt
    using Statistics
    include("modelling_funcs.jl")
    export generate_dynamics, dynamics_acc_ctrl, pend_sim,dynamics_acc_ctrl_terms

    include("plotting_funcs.jl")
    export plot_FFT, rot_pendulum_animator, plot_rot_pendulum, general_lin_interp

    include("TrajOptim.jl")
    export  SpinUpTrajectory
end # module RotaryInvertedPendulum
