__precompile__()
module RotaryInvertedPendulum

    using Symbolics
    using LinearAlgebra
    using DifferentialEquations
    using Plots
    using Interpolations

    include("TrajOptim.jl")
    export  SpinUpTrajectory
end # module RotaryInvertedPendulum
