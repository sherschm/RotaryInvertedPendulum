# Rotary Pendulum Simulation

- This repo contains a Julia script that derives the nonlinear state-space equations of motion of the EduKit inverted rotary pendulum,
intended to be used in the UKAEA PACE training programme. The derivation of this model is detailed in [the pdf.](https://github.com/sherschm/RotaryInvertedPendulum/blob/main/Modelling%20%26%20Simulation%20of%20a%20rotary%20inverted%20pendulum.pdf)

- The derived model can be simulated and animated:

![uncontrolled cartpole gif](./rotary_pendulum_anim.gif)

![pendulum response](./response.png)

## Preliminaries
[Install Julia](https://docs.julialang.org/en/v1/manual/installation/)

From a command prompt, run Julia

```bash
julia
```
In the Julia REPL, import the required packages by running:
```bash
using Pkg
Pkg.add("Symbolics")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("DifferentialEquations")
Pkg.add("Interpolations")
exit()
```

## To run the modelling & simulation script from command line:
Clone the repository and move to its directory.

Run the modelling script:
```bash
julia main.jl
```
This commented script runs through the model derivation and simulation.

## Next steps...
Control!
