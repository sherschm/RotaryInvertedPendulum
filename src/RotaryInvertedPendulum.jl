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
    
    include("TrajOptim.jl")
    export  SpinUpTrajectory
end # module RotaryInvertedPendulum
