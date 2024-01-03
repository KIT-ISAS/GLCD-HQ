push!(LOAD_PATH, pwd())

# Dependencies:
import Pkg; Pkg.add.(["Revise","QuadGK","ForwardDiff","Optimization","OptimizationOptimJL","CairoMakie"])

using Revise
using LCDGauss
import OptimizationOptimJL
using CairoMakie; CairoMakie.activate!()


# Gaussian Parameters
L = 30       # number of samples
σ = [0.5,1]  # standard deviations 

# Solver Parameters
bmin = 0.1; bmax = 10 # Steinbring C++, for standard-normal
# bmin = 0.0001; bmax = 100 # UDH Matlab
solver = OptimizationOptimJL.ConjugateGradient() # OptimizationOptimJL.ConjugateGradient(), OptimizationOptimJL.BFGS()
solver_reltol = 1e-6
quad_reltol = 1e-9
fast = true

# Compute Samples
x = LCDGauss.sample(σ; L, bmin, bmax, solver, solver_reltol, quad_reltol, fast) 

# Plot
fig = Figure(size = (700,700), figure_padding = 1)
ax = Axis(fig[1, 1], aspect=DataAspect())
h = scatter!(ax, x[1,:], x[2,:]) 
h.markersize[] = h.markersize[]*3

display(fig); # opens PNG viewer

save("gaussian.pdf", fig, pt_per_unit = 1) # saves PDF

