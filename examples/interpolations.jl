using Pkg
Pkg.activate("examples/")

using TransitionFunctionIteration
using StructArrays
using Plots 

# Define a 2D grid of points
ssg = StateSpaceGrid(x=-5.:5., y=-3.:3.)

# Calculate two functions to interpolate
Z = [pt.x^2 + pt.y^3 for pt in points(ssg)]
W = [pt.x - pt.y for pt in points(ssg)]

# Create interpolant
fct = Interpolant(ssg, StructArray(( ; Z, W)))

# Plot 2D
p1 = plot(fct, :W, (;))

# Plot 1D
p2 = plot(fct, :Z, (; y=-1.0))

# Lay out side by side 
plot(p1, p2, layout=(1,2), size=(800,400))