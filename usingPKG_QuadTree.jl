# Documentation for module:
# https://docs.juliahub.com/Quadtrees/dYG06/0.1.0/

import Pkg; Pkg.add("Quadtrees")
using Quadtrees

mutable struct Nodes
    x::Real
    y::Real
    radius::Real
end

Quadtrees.position(p::Nodes) = (p.x, p.y)

qt = Quadtree{Nodes}((0, 0), 400, 400, 4)

for i in 0:999
    # two random numbers between 0 and 200
    x, y = rand(Float64, 2) * 200
    # add 1 to make sure the radius is not zero
    radius = rand(Float64) * 49 + 1
    p = Particle(x, y, radius)
    # insert the particle into the quadtree.
    insert!(qt, p)
end