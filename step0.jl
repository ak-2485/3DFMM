
"""
"""

using GridStruct
using BoxStruct
using ListConstructors
using OctreeConstructor
using CoefHelpers

function step0(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
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

    Calculate the sum of all local contributions in each leaf box.
    """
    ψ = 0
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        if !haschildren(box)
            ψ += evallocalpotential(box, grid, p)
        end
    end
    return ψ
end # step 6

function step7(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 7

    Calculate near neighbor contributions directly.
    """
    ψ = 0
    boxids = keys(grid.boxes)
    for boxid in boxids
        box = grid.boxes[boxid]
        particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()
        if !haschildren(box)
            L1boxids = box.L1
            for L1boxid in L1boxids
                L1box = grid.boxes[L1boxid]
                L1boxparticledict = getparticledict(L1box, grid)
                particledict = merge(particledict,L1boxparticledict)
            end
            ψ += evaldirectpotential(box, grid, p, particledict)
        end
    end

    return ψ
end # step 7

function step8(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 8

    Calculate nearest neighbor contributions directly.
    """
    ϕ = 0
    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
        particledict = Dict()
        if !haschildren(box)
            L3boxids = box.L3
            for L3boxid in L3boxids
                L3box = grid.boxes[L3boxid]
                L3boxparticledict = getparticledict(L3box, grid)
                if length(L3boxparticledict) > p^2
                    ϕ += evalmultipolepotential(box,L3box,grid,p)
                else
                    ϕ += evaldirectpotential(box, grid, p, L3boxparticledict)
                end
            end
        end
    end

    return ϕ
end # step 8

function FMM3D(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
    p::Int64, nmin::Int64)
     """
     """
     grid = step0(particles,nmin)
     step1(grid,p)
     step2(grid,p)
     step3(grid,p)
     step4(grid,p)
     step5(grid,p)
     ψ1 = step6(grid,p)
     ϕ2 = step7(grid,p)
     ϕ3 = step8(grid,p)

     Φ = ψ1 + ϕ2 + ϕ3
     println(Φ)
end

particles = Dict([(1,((0.24,0.22,0.26),3.0)),(2,((0.72,0.76,0.72),4.0)),
(3,((0.54,0.25,0.25),3.0)),(4,((0.73,0.65,0.72),2.0)), (5,((0.12,0.22,0.22),2.0)),
(6,((0.72,0.72,0.14),3.0))])
FMM3D(particles,4,2)
