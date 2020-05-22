"""
"""
module EnsembleTests

using Random, Distributions

export particledist2, particledist1, particledist3, rndmparticledist3
export particledist1ex, clusterdist

function particledist1()
    """
    """
    particles = Dict([(1,((0.24,0.22,0.26),3.0)),(2,((0.72,0.76,0.72),4.0)),
    (3,((0.54,0.25,0.25),3.0)),(4,((0.73,0.65,0.72),2.0)), (5,((0.12,0.22,0.22),2.0)),
    (6,((0.72,0.72,0.14),3.0))])

    minbound = [0.0,0.0,0.0]
    maxbound = [1.0,1.0,1.0]

    return particles, minbound, maxbound
end

function clusterdist()
    """
    """
    particles = Dict([(1,((-4.5,-4.5,-4.5),3.0)),(2,((4.5,-4.5,-4.5),4.0)),
    (3,((-4.5,4.5,-4.5),3.0)),(4,((-4.5,-4.5,4.5),2.0)), (5,((4.5,4.5,4.5),2.0)),
    (6,((4.5,4.5,-4.5),3.0))])

    minbound = [-10.0,-10.0,-10.0]
    maxbound = [10.0,10.0,10.0]

    return particles, minbound, maxbound
end

function particledist1ex()
    """
    """
    coordsx = [x for x in -9.5:3:9.5 for y in -9.5:3:9.5 for z in -9.5:3:9.5]
    coordsy = [y for x in -9.5:3:9.5 for y in -9.5:3:9.5 for z in -9.5:3:9.5]
    coordsz = [z for x in -9.5:3:9.5 for y in -9.5:3:9.5 for z in -9.5:3:9.5]

    minbound = [-10.0,-10.0,-10.0]
    maxbound = [10.0,10.0,10.0]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()

    sz=length(coordsx)
    coords = zeros(sz,3)
    coords[:,1] = coordsx
    coords[:,2] = coordsy
    coords[:,3] = coordsz

    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = 0.25*(-1)^i *(x+y)*z
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, particledict, minbound, maxbound
end


function particledist2()
    """
    """
    coordsx = [x for x in -9.875:0.5:9.875 for y in -9.875:0.5:9.875 for z in -9.875:0.5:9.875]
    coordsy = [y for x in -9.875:0.5:9.875 for y in -9.875:0.5:9.875 for z in -9.875:0.5:9.875]
    coordsz = [z for x in -9.875:0.5:9.875 for y in -9.875:0.5:9.875 for z in -9.875:0.5:9.875]

    minbound = [-10.0,-10.0,-10.0]
    maxbound = [10.0,10.0,10.0]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()

    sz=length(coordsx)
    coords = zeros(sz,3)
    coords[:,1] = coordsx
    coords[:,2] = coordsy
    coords[:,3] = coordsz

    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = rand(-5:0.1:5)
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, particledict, minbound, maxbound
end

function particledist3()
    """
    """
    coordsx = [x for x in -0.875:0.5:0.875 for y in -0.875:0.5:0.875 for z in -0.875:0.5:0.875]
    coordsy = [y for x in -0.875:0.5:0.875 for y in -0.875:0.5:0.875 for z in -0.875:0.5:0.875]
    coordsz = [z for x in -0.875:0.5:0.875 for y in -0.875:0.5:0.875 for z in -0.875:0.5:0.875]

    minbound = [-1.0,-1.0,-1.0]
    maxbound = [1.0,1.0,1.0]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()

    sz=length(coordsx)
    coords = zeros(sz,3)
    coords[:,1] = coordsx
    coords[:,2] = coordsy
    coords[:,3] = coordsz

    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = rand(-2:0.2:2)
        #q = (-1)^i*(x+y)/z
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, particledict, minbound, maxbound
end

function rndmparticledist3(num::Int64)
    """
    """
    minbound = [-0.5,-0.5,-0.5]
    maxbound = [0.5,0.5,0.5]

    coords = rand(Uniform(-0.5,0.5), num,3)
    sz=size(coords)[1]

    particledict = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64,Float64}}()
    particles = Dict{Int64,Tuple{Tuple{Float64,Float64,Float64},Float64}}()
    for i = 1:sz
        x,y,z = coords[i,:]
        coord = (x,y,z)
        q = rand(-5:0.1:5)
        particledict[i] = (coord,q,0.0)
        particles[i] = (coord,q)
    end

    return particles, minbound, maxbound
end

end #module
