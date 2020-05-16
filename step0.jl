
"""
"""

using GridStruct
using BoxStruct
using ListConstructors
using OctreeConstructor
using CoefHelpers

function step0(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}},
    nmin::Int64)
    """
    Build Tree Structure and create Lists

    input: particles as dictionary

    preconditions:

    returns: an adaptive grid for computing the potential due to all particles
    in the dictionary.

    """
    grid = Grid()
    grid.nmin = nmin
    grid.maxbound = [1,1,1]
    grid.minbound = [0,0,0]
    grid.particles = particles
    grid.center = [1,1,1] .* (grid.maxbound .+ grid.minbound)[1]/2
    grid.size = grid.maxbound[1] .- grid.minbound[1]
    grid.numparticles = length(grid.particles)
    octree(grid)
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        L1!(box, grid.boxes)
        L2!(box, grid.boxes)
        L3!(box, grid.boxes)
    end
    for boxid in boxids
        box = grid.boxes[boxid]
        L4!(box, grid.boxes)
    end

    return grid
end

function step1(grid::Grid, p::Int64)
    """
    Upward Pass, Step 1

    Calculate the multipole coefficients for all leaf boxes.
    """
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        if !haschildren(box)
            mcoef!(box,p,grid)
        end
    end
end

function step2(grid::Grid, p::Int64)
    """
    Upward Pass, Step 2

    Starting at the finest level, translate the multipole coefficients of
    children to parents.
    """
    maxlevel = length(grid.levels)-1
    for i = maxlevel:-1:2
        boxesonlevel = grid.levels[i]
        for boxid in boxesonlevel
            child = grid.boxes[boxid]
            Mjk = mcoeftrans(child,p)
            parent = grid.boxes[child.parent]
            parent.multipole_coef = Mjk
        end #loop over boxes on level
    end #loop over levels
end #step 2

function step3(grid::Grid, p::Int64)
    """
    Downward Pass, Step 3

    For each box b, find the local expansion coefficients due to all charges
    in L4(b).
    """
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        L4particles = Set{Int64}()
        # Collect the particle ids of the particles in all boxes in L4 of box
        for L4boxid in box.L4
            L4box = grid.boxes[L4boxid]
            union!(L4particles,L4box.particles)
        end
            lcoef!(box,p,L4particles,grid)
    end #boxes
end #step 3

function step4(grid::Grid, p::Int64)
    """
    Downward Pass, Step 4

    For each box b, convert the multipole expansions of all boxes in L2(b)
    to local expansions about the center of b and add the result to ψ(b).
    """
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        for boxid in box.L2
            L2box = grid.boxes[boxid]
            Ljk = mtolconversion(L2box, p)
            lcoefsum!(box, Ljk, p)
        end
    end
end #step 4

function step5(grid::Grid, p::Int64)
    """
    Downward Pass, Step 5

    For level l, for each box b on l, translate the local expansion of b
    to its children.
    """
    levels = sort(collect(keys(grid.levels)))
    for levelid in levels
        boxesonlevel = grid.levels[levelid]
        for boxid in boxesonlevel
            box = grid.boxes[boxid]
            for boxid in box.children
                child = grid.boxes[boxid]
                Ljk = lcoeftrans(box, p)
                lcoefsum!(child, Ljk, p)
            end
        end
    end
end # step 5

function step6(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 6

    Calculate local contributions on finest level for each particle.
    """

    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
        if !haschildren(box)
            localatpoints!(box, grid, p)
        end
    end

end # step 6

function step7(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 7

    Calculate near neighbor contributions on finest level at each particle directly.
    """

    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
        if !haschildren(box)
            for L1boxid in box.L1
                L1box = grid.boxes[L1boxid]
                directatpoints!(box, L1box, grid, p)
            end
        end
    end

end # step 7

function step8(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 8

    Calculate semi-near neighbor multipole or direct contributions on finest level
    at each particle.
    """
    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
        if !haschildren(box)
            for L3boxid in box.L3
                L3box = grid.boxes[L3boxid]
                if length(L3box.particles) > p^2
                    multipoleatpoints!(box,L3box,grid,p)
                else
                    directatpoints!(box,L3box,grid,p)
                end
            end
        end
    end

end # step 8

function FMM3D(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
    p::Int64, nmin::Int64)
     """
     """
     particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}()
     for (key,value) in particles
         coords = value[1]
         q = value[2]
         particledict[key] = (coords,q,0.0)
     end
     grid = step0(particledict,nmin)
     step1(grid,p)
     step2(grid,p)
     step3(grid,p)
     step4(grid,p)
     step5(grid,p)
     step6(grid,p)
     step7(grid,p)
     step8(grid,p)
end

particles = Dict([(1,((0.24,0.22,0.26),3.0)),(2,((0.72,0.76,0.72),4.0)),
(3,((0.54,0.25,0.25),3.0)),(4,((0.73,0.65,0.72),2.0)), (5,((0.12,0.22,0.22),2.0)),
(6,((0.72,0.72,0.14),3.0))])
FMM3D(particles,4,2)
