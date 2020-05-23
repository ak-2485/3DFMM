"""
"""
module FMM3D

using FMMSteps

export fastmultipole

function fastmultipole(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
    minbound::Array{Float64},maxbound::Array{Float64},nmax::Int64, p::Int64)
     """
     """
     grid = step0(particles,minbound,maxbound,nmax,p)
     println("maximum number of particles in a leaf box:", nmax)
     println("expansion order:", p)
     println("number of particles:", length(particles))
     println("number of levels: ", length(grid.levels))
     println("number of boxes: ", length(grid.boxes))
     step1(grid,p)
     step2(grid,p)
     step3(grid,p)
     step4(grid,p)
     step5(grid,p)
     step6(grid,p)
     step7(grid,p)
     step8(grid,p)

     return grid.particles
end

end #module
