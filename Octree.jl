
"""

"""
module OctreeConstructor

using treehelpers

export Octree

function Octree(points::Dict{Int64,Tuple{Float64,Float64,Float64}},
        bounds::Array{Array{Float64,1},1},
        nmin::Int)
        """

        input:
            points: a dictinoary of particles in the domain with particle
                index, numbered (1...total-number-of-particles),
                 as key and particle location as value.
            bounds: bounds of the computational domain, the first array index
                is the least vertex point of the domain, the second array index
                is the greatest vertex point of the domain.
            nmin : boxes with fewer than nmin particles won't be refined.
        preconditions:
            bounds: a non-empty array of 2 non-empty arrays
            nmin: must be greater than 1 for recursion to halt 
        returns:
            levelDict : A dictionary with level number as key and box index as value
            treeDict : A dictionary with box index as key and box as value
        """

    # A stack to keep track of the total number of box indices
    indexCount = Stack{Int}()
    # A stack to keep track of the total number of levels
    levelCount = Stack{Int}()
    # A dictionary for keeping track of recursive level visitation
    visited = Dict{Int64,Bool}()
    # A dictionary with box index as key and box as value
    treeDict = Dict{Int64,Box}()
    # A dictionary with level number as key and box index as value
    levelDict = Dict{Int64,Array{Int64,1}}()
    # An array for calculating vertices and centers
    # arranged s.t. |7,8|5,6|
    #              |1,2|3,4|
    s = [[-1,-1,1],[-1,-1,-1],[1,-1,1],[1,-1,-1],
            [1,1,1],[1,1,-1],[-1,1,1],[-1,1,-1]]

    l = length(points)
    # Create the initial Box, which is the entire grid
    B = Box()
    B.num_particles = l
    B.particle_list = collect(1:l)
    B.center = [1,1,1] .* (bounds[2] .+ bounds[1])[1]/2
    B.size = bounds[2][1] .- bounds[1][1]
    for i in 1:8
        B.vertices[i] = 0.5 * B.size .* s[i] + B.center
    end
    B.max_bound = B.vertices[5]
    B.min_bound = B.vertices[2]
    # Begin populating the tree and level dictionaries
    treeDict[B.index] = B
    levelDict[B.level] = [B.index]
    # Only divide if sufficient points are available
    if l < nmin return end
    # Recursively divide the box
    OctreeDivide(B,points,nmin)
    # update colleagues
    levels = keys(levelDict)
    for level in levels
        getColleagues(treeDict, levelDict, level)
    end
    return levelDict, treeDict
end

function OctreeDivide(Parent::Box,
    points::Dict{Int64,Tuple{Float64,Float64,Float64}}, nmin::Int64)
    """

    input:
    returns : none
    """
    # If the refinement condition has been met, exit
    if Parent.num_particles < nmin return end
    # If the level hasn't been visited, visit it
    if !haskey(visited,Parent.level)
        push!(levelCount,1)
        visited[Parent.level] = true
        levelDict[length(levelCount)] = []
    end
    # Recursively make babies
    Children = [Box() for i in 1:8]
    j = 1
    for child in Children
        child.size = 0.5 * Parent.size
        child.center = Parent.center + 0.5 * child.size * s[j]
        for i in 1:8
            child.vertices[i] = 0.5 * child.size .* s[i] + child.center
        end
        child.max_bound = child.vertices[5]
        child.min_bound = child.vertices[2]
        child.particle_list = pointsIn(points,child)
        child.num_particles = length(child.particle_list)
        # Ignore empty children
        if child.num_particles > 0
            child.parent = Parent.index
            child.index = length(indexCount) + 1
            push!(indexCount,1)
            treeDict[child.index] = child
            push!(Parent.children,child.index)
            child.level = Parent.level + 1
            levelDict[child.level] = [levelDict[child.level]; child.index]
            OctreeDivide(child,points,nmin)
        end
        j = j + 1
    end
end
