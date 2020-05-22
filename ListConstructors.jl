
"""
"""
module ListConstructors

using BoxStruct

export L1!, L2!, L3!, L4!, dirlists!

function L1!(box1::Box, boxdict::Dict{Int64,Box})
    """
    Populates box1.L1: The list of box1's leaf colleagues. If box1 is a parent
    box, then box1.L1 is empty.

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: none
    """
    L1 = Set{Int64}()
    # If box1 is a parent box, l1 is empty
    if haschildren(box1) return end
    boxes = keys(boxdict)
    for bxid in boxes
        box2 = boxdict[bxid]
        if !haschildren(box2) && areadjacent(box1,box2)
            push!(L1,box2.index)
        end
    end
    box1.L1 = L1
end


function L2!(box1::Box, boxdict::Dict{Int64,Box})
    """
    Populates box1.L2: a set of all of the children of the colleagues of
    box1's parent that are well separated from box1.

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: none
    """
    # Get the colleauges of box1's Parent
    Parent = boxdict[box1.parent]
    colParent = Parent.colleagues
    allchildren = Set{Int64}()
    L2 = Set{Int64}()
    # Get a set of all children of the colleagues of box1's parent
    for bxid in collect(colParent)
        box2 = boxdict[bxid]
        allchildren = union(allchildren, box2.children)
    end
    # Collect the boxids from allchildren that are well separated from box1
    for i in allchildren
        box3 = boxdict[i]
        if arewellseparated(box1,box3)
            push!(L2,box3.index)
        end
    end
    box1.L2 = L2
end

function dirlists!(box1::Box, boxdict::Dict{Int64,Box})
    """
    A subset of box1.L2 of boxes that are above box1 (+z dir) and are
    separated from box1.L2 by at least one box.

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: none
    """
    uplist = Set{Int64}()
    downlist = Set{Int64}()
    northlist = Set{Int64}()
    southlist = Set{Int64}()
    eastlist = Set{Int64}()
    westlist = Set{Int64}()
    # Get the maximum z value of box1's vertices
    zmaxb1 = box1.max_bound[3]
    zminb1 = box1.min_bound[3]
    ymaxb1 = box1.max_bound[2]
    yminb1 = box1.min_bound[2]
    xmaxb1 = box1.max_bound[1]
    xminb1 = box1.min_bound[1]
    for bxid in box1.L2
        box2 = boxdict[bxid]
        # Get the minimum z value of vertices for a box in box1.L2
        zmaxb2 = box2.max_bound[3]
        zminb2 = box2.min_bound[3]
        ymaxb2 = box2.max_bound[2]
        yminb2 = box2.min_bound[2]
        xmaxb2 = box2.max_bound[1]
        xminb2 = box2.min_bound[1]
        if arewellseparated(box1,box2)
            if zminb2 > zmaxb1
                push!(uplist,bxid)
            end
            if zmaxb2 < zminb1
                push!(downlist,bxid)
            end
            if yminb2 > ymaxb1
                push!(northlist,bxid)
            end
            if ymaxb2 < yminb1
                push!(southlist, bxid)
            end
            if  xminb2 > xmaxb1
                push!(eastlist, bxid)
            end
            if xmaxb2 < xminb1
                push!(westlist, bxid)
            end
        end
    end
    box1.uplist = uplist
    box1.downlist = downlist
    box1.northlist = northlist
    box1.southlist = southlist
    box1.eastlist = eastlist
    box1.westlist = westlist
end

function getdescendants(box1::Box, boxdict::Dict{Int64,Box})
    """
    Recursively find all descendants of a box1

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: a set of boxids of the descendants of box1.
    """
    visited = Set{Int64}()
    function descend(box1::Box)
        # Visit
        push!(visited, box1.index)
        # Base
        if !haschildren(box1) return visited end
        # Recurse
        children = collect(box1.children)
        for child in children
            box2 = boxdict[child]
            descend(box2)
        end
        return visited
    end
    return descend(box1)
end #getdescendants

function L3!(box1::Box, boxdict::Dict{Int64,Box})
    """
    Populates box1.L3: The list of all descendents of box1's colleagues that are
    not adjacent to box1. If box1 is a parent box, L3 is empty. Any box in L3
    is smaller than box1.

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: none.
    """
    L3 = Set{Int64}()
    # If box1 is a parent box, L3 is empty
    if haschildren(box1) return L3 end
    # get all descendants of all colleagues of box1
    for colid in collect(box1.colleagues)
        cBox = boxdict[colid]
        L3 = union(L3,getdescendants(cBox,boxdict))
    end
    for boxid in collect(L3)
        box2 = boxdict[boxid]
        # remove those that are adjacent to box1
        if areadjacent(box1,box2)
            pop!(L3, boxid)
        end
        # remove those whose parents aren't adjacent to box1
        box3 = boxdict[box2.parent]
        if !areadjacent(box2,box3)
            pop!(L3, box3.index)
        end
    end
    box1.L3 = L3
end # L3!

function L4!(box1::Box, boxdict::Dict{Int64,Box})
    """
    The list of all boxes c such that c∈L3(box1). All boxes in L4(box1) are
    childless and larger than b.

    inputs:
        box1: the box under consideration.
        boxdict: a dictionary of boxes from the grid.
    preconditions: box1 must be an element of boxdict; boxdict must be the
        entire dictionary of boxes in a grid.
    returns: none.
    """
    L4 = Set{Int64}()
    boxids = keys(boxdict)
    for boxid in boxids
        box2 = boxdict[boxid]
        if box1.index in box2.L3
            push!(L4,boxid)
            @assert isempty(box2.children) ["boxes added to box.L4 are childless"]
            @assert box2.level<box1.level ["boxes added to box.L4 are larger than box"]
        end
    end
    box1.L4 = L4
end # L4!

#end module
end
