"""
Due to compile time for Plots, plotting for test ensembles is kept separate
from the test themselves.
"""

using GridStruct
using BoxStruct
using Plots
using EnsembleTests

function plotgrid(grid::Grid)
    xs = zeros(8,1)
    ys = zeros(8,1)
    zs = zeros(8,1)
    flag = false
    plt = 0
    for (boxid,val) in grid.boxes
        box = grid.boxes[boxid]
        if !(BoxStruct.haschildren(box))
            v = box.vertices
            for i = 1:8
                xs[i] = v[i][1]
                ys[i] = v[i][2]
                zs[i] = v[i][3]
            end
            plt = plot!([xs[2],xs[1]],[ys[2],ys[1]],[zs[2],zs[1]],legend = false, axis=nothing, foreground_color_subplot=colorant"white")
            plot!([xs[2],xs[8]],[ys[2],ys[8]],[zs[2],zs[8]])
            plot!([xs[2],xs[4]],[ys[2],ys[4]],[zs[2],zs[4]])

            plot!([xs[6],xs[8]],[ys[6],ys[8]],[zs[6],zs[8]])
            plot!([xs[6],xs[5]],[ys[6],ys[5]],[zs[6],zs[5]])
            plot!([xs[6],xs[4]],[ys[6],ys[4]],[zs[6],zs[4]])

            plot!([xs[7],xs[8]],[ys[7],ys[8]],[zs[7],zs[8]])
            plot!([xs[7],xs[5]],[ys[7],ys[5]],[zs[7],zs[5]])
            plot!([xs[7],xs[1]],[ys[7],ys[1]],[zs[7],zs[1]])

            plot!([xs[3],xs[5]],[ys[3],ys[5]],[zs[3],zs[5]])
            plot!([xs[3],xs[4]],[ys[3],ys[4]],[zs[3],zs[4]])
            plot!([xs[3],xs[1]],[ys[3],ys[1]],[zs[3],zs[1]])

            flag = true
        end
    end
    coordsx = zeros(length(grid.particles))
    coordsy = zeros(length(grid.particles))
    coordsz = zeros(length(grid.particles))
    i = 1
    for (partid,vals) in grid.particles
        (x,y,z) = vals[1]
        coordsx[i] = x
        coordsy[i] = y
        coordsz[i] = z
        i += 1
    end
    scatter!(coordsx, coordsy, coordsz, axis=nothing)
    if flag == true display(plt) end

end #function

particles, _ ,minbound,maxbound = particledist2()
nmax = 10
p = 4
grid = step0(particles,minbound,maxbound,nmax)
plotgrid(grid)
