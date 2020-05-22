
"""
"""
module FMMSteps

using GridStruct
using BoxStruct
using ListConstructors
using OctreeConstructor
using CoefHelpers

export step0, step1, step2, step3, step4, step5, step6, step7, step8

function step0(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}},
    minbound::Array{Float64,1}, maxbound::Array{Float64,1},nmax::Int64,p::Int64)
    """
    Build Tree Structure and create Lists

    input: particles as dictionary

    preconditions:

    returns: an adaptive grid for computing the potential due to all particles
    in the dictionary.

    """
    grid = Grid()
    grid.nmax = nmax
    grid.maxbound = maxbound
    grid.minbound = minbound

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}}()
    for (key,value) in particles
        coords = value[1]
        q = value[2]
        particledict[key] = (coords,q,0.0)
    end

    grid.particles = particledict
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
    lmax = (p+1)^2
    for boxid in boxids
        box = grid.boxes[boxid]
        L4!(box, grid.boxes)
        box.multipole_coef = zeros(ComplexF64,lmax)
        box.local_coef = zeros(ComplexF64,lmax)
    end
    println("level dict:", grid.levels)

    return grid
end

function step1(grid::Grid, p::Int64)
    """
    Upward Pass, Step 1

    Calculate the multipole coefficients for all leaf boxes.
    """
    for boxid in keys(grid.boxes)
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
    for i = maxlevel:-1:1
        boxesonlevel = grid.levels[i]
        for boxid in boxesonlevel
            child = grid.boxes[boxid]
            parent = grid.boxes[child.parent]
            Mjk = mcoeftrans(child,parent,p)
            parent.multipole_coef .+= Mjk
        end #loop over boxes on level
    end #loop over levels
end #step 2

function step3(grid::Grid, p::Int64)
    """
    Downward Pass, Step 3

    For each box b, find the local expansion coefficients due to all charges
    in L4(b).
    """
    for boxid in keys(grid.boxes)
        box = grid.boxes[boxid]
        for l4boxid in box.L4
            l4box = grid.boxes[l4boxid]
            if length(box.particles) > p^2
                println("num particles in L4s > p^2")
                lcoef!(box,l4box,p,grid)
            else
                directatpoints!(box, l4box, grid, p)
            end
         end
    end

end #step 3

function step4(grid::Grid, p::Int64)
    """
    Downward Pass, Step 4

    For level l, for each box b on l, convert the multipole expansions of all
    boxes in L2(b) to local expansions about the center of b and add the result to b's local
    expansion.
    """
    levels = sort(collect(keys(grid.levels)))
    for levelid in levels
        boxesonlevel = grid.levels[levelid]
        for boxid1 in boxesonlevel
            box = grid.boxes[boxid1]
            for boxid2 in box.L2
                L2box = grid.boxes[boxid2]
                Ljk = mtolconversion(L2box, box, p)
                box.local_coef .+= Ljk
            end
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
        for boxid1 in boxesonlevel
            box = grid.boxes[boxid1]
            for boxid2 in box.children
                child = grid.boxes[boxid2]
                Ljk = lcoeftrans(box, child, p)
                child.local_coef .+= Ljk
            end
        end
    end
end # step 5

function step6(grid::Grid, p::Int64)
    """
    Evaluate Potentials, Step 6

    Calculate local contributions on finest level at each particle.
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
                    println("multipoleatpoints!")
                    multipoleatpoints!(box,L3box,grid,p)
                else
                    directatpoints!(box,L3box,grid,p)
                end
            end
        end
    end

end # step 8

end #module
