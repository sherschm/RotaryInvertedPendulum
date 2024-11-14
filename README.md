# Julia Rotary Pendulum Simulation

- Derives the nonlinear equations of motion of the EduKit inverted rotary pendulum,
intended to be used in the PACE training programme.

- Simulates and animates the equations of motion:

![uncontrolled cartpole gif](./rotary_pendulum_anim.gif)

![pendulum response](./response.png)

For the details on the derivation, take a look at the pdf

## Preliminaries
[Install Julia](https://docs.julialang.org/en/v1/manual/installation/)

In Julia, import the required packages by running:
```bash
using Pkg
Pkg.add("Symbolics")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("DifferentialEquations")
Pkg.add("Interpolations")
```

## Installation:

From a command prompt, navigate to your chosen clone folder, then run

```bash
git clone git@git.ccfe.ac.uk:sherschm/rotaryinvertedpendulum.git
```
```bash
cd rotaryinvertedpendulum
```
## To run the modelling & simulation script from cmd prompt:
```bash
julia main.jl
```


