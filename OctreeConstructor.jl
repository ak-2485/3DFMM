
"""

"""
module OctreeConstructor

using DataStructures
using GridStruct
using BoxStruct
using TreeHelpers

export octree

function octree(Grid::Grid)
        """
        Build Tree Structure by refining the domain adaptively on boxes that contain
        more than nmax particles.
        input:
            points: a dictinoary of particles in the domain with particle
                index, numbered (1...total-number-of-particles),
                 as key and particle location and strength as value.
            bounds: bounds of the computational domain, the first array index
                is the least vertex point of the domain, the second array index
                is the greatest vertex point of the domain.
            nmax : boxes with fewer than nmax particles won't be refined.
        preconditions:
            bounds: a non-empty array of 2 non-empty arrays
            nmax: must be greater than 1 for recursion to halt
        returns:
            leveldict : A dictionary with level number as key and box index as value
            boxdict : A dictionary with box index as key and box as value
    """
    # keep track of the total number of box indices in order to set
    # new box indices
    indexacc = 0
    # keep track of the total number of levels in order to set
    # new box levels
    levelacc = 0
    # Keep track of levels visited
    visited = Set{Int64}()
    s = Grid.signarray
    # Create the initial Box, which is the entire grid
    B = Box()
    B.numparticles = Grid.numparticles
    B.particles = keys(Grid.particles)
    B.center = Grid.center
    B.size = Grid.size
    for i in 1:8
        B.vertices[i] = 0.5 * B.size .* Grid.signarray[i] + B.center
    end
    B.max_bound = Grid.maxbound
    B.min_bound = Grid.minbound
    # Begin populating the tree and level dictionaries
    Grid.boxes[0] = B
    Grid.levels[0] = Set(B.index)
    # Only divide if sufficient points are available
    if B.numparticles <= Grid.nmax return end
    # Recursively divide the box is it contains more than nmax particles
    function octreedivide(Parent::Box)
        """

        input:
        returns : none
        """
        # If the refinement condition has been met, exit
        if Parent.numparticles <= Grid.nmax return end
        # If the level hasn't been visited, visit it
        if !(Parent.level in visited)
            levelacc += 1
            push!(visited,Parent.level)
            Grid.levels[levelacc] = Set()
        end
        # Recursively make babies
        Children = [Box() for i in 1:8]
        j = 1
        # Create set to use for allocation of particles; to ensure particles
        # is located in only one box
        particledict = copy(Grid.particles)
        for child in Children
            child.size = 0.5 * Parent.size
            child.center = Parent.center + 0.5 * child.size * s[j]
            for i in 1:8
                child.vertices[i] = 0.5 * child.size .* s[i] + child.center
            end
            child.max_bound = child.vertices[5]
            child.min_bound = child.vertices[2]
            child.particles = particlesin(particledict,child)
            for particle in child.particles
                pop!(particledict,particle)
            end
            child.numparticles = length(child.particles)
            # Ignore empty children
            if child.numparticles > 0
                child.parent = Parent.index
                indexacc += 1
                child.index = indexacc
                Grid.boxes[child.index] = child
                push!(Parent.children,child.index)
                child.level = Parent.level + 1
                if haskey(Grid.levels,child.level)
                    Grid.levels[child.level] = push!(Grid.levels[child.level],child.index)
                else
                    Grid.levels[child.level] = Set(child.index)
                end
                octreedivide(child)
            end
            j = j + 1
        end
        # update colleagues
        boxids = keys(Grid.boxes)
        for boxid in boxids
            colleagues!(Grid.boxes[boxid],Grid.levels,Grid.boxes)
        end
    end
    return octreedivide(B)
end

end #module OctreeConstructor
