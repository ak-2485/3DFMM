
"""
"""
module TreeHelpers

using BoxStruct

export colleagues!, domainBounds

function colleagues!(box1::Box, leveldict::Dict{Int64,Set{Int64}},
    boxdict::Dict{Int64,Box})
    """

    """
    level = box1.level
    # set of box ids on box1's level
    boxids = leveldict[level]
    for boxid in boxids
        box2 = boxdict[boxid]
        if iscolleague(box1,box2)
            push!(box1.colleagues,box2.index)
            push!(box2.colleagues,box1.index)
        end
    end
end

function domainBounds(particles::Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Complex{Float64}}})
    """
    This function can be used to determine a box containing all particles in
    a given particle dictionary.

    inputs: particles as a dictionary

    returns: minimum and maximum boundary of a box that contains all particles

    """
    sortpartArr = sort(collect(particles),by=x->x[2])
    len = length(particles)
    return sortpartArr[1][2][1], sortpartArr[len][2][1]
end

end #module
