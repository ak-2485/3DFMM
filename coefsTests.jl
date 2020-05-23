
using Test
using LinearAlgebra
using CoefHelpers
using GridStruct
using EnsembleTests
using BoxStruct
#using PlotGrid
import FMMSteps
import FMM3D

function mtomTest4()
    """
    Test the M2M operator by first comparing the moments of a parent
    box to the moments of the same parent box after translation from
    children.
    """
    println("Testing m2m: parent <--- children")
    particles, minbound, maxbound = particledist1()
    p = 10
    println("p = ", p)
    nmax = length(particles)
    grid1 = FMMSteps.step0(particles,minbound,maxbound,nmax,p)
    # set the moments of the parent box when it is childless
    FMM3D.step1(grid1,p)
    # on a new grid, set the moments of the children then pass them to parent
    nmax = 2
    grid2 = FMMSteps.step0(particles,minbound,maxbound,nmax,p)
    FMM3D.step1(grid2,p)
    FMM3D.step2(grid2,p)
    # compare moments of the parent
    parent1 = grid1.boxes[0]
    parent2 = grid2.boxes[0]
    Mparent1 = parent1.multipole_coef
    Mparent2 = parent2.multipole_coef
    num = 0
    den = 0
    for i = 1:length(Mparent1)
        Mnmexpected = Mparent1[i]
        Mnm = Mparent2[i]
        num += abs(Mnmexpected - Mnm)^2
        den += abs(Mnmexpected)^2
    end
    E = sqrt(num/den)
    @test E < 10^-6

end

mtomTest4()

function ltolTest1()
    """
    Test the local to L2L operator by comparing the lcoefs of the parent due to
    a colleague translated to a child and the lcoefs of the child due to the
    same colleague of the parent. These values are equal for p=0.
    """
    println("Testing l2l: parent --> child due to parent colleague")
    particles,_,minbound, maxbound = particledist1ex()
    # This test requires p = 0 !
    p = 0
    nmax = 20
    grid = FMM3D.step0(particles,minbound,maxbound,nmax,p)
    # pick a parent box from the grid
    parent = grid.boxes[1]
    # get a colleague
    colparent = grid.boxes[19]
    @assert parent.level == colparent.level
    # set the local coefs of the parent due to colleague's particles
    lcoef!(parent,colparent,p,grid)
    # get array of coefs that would pass to the child
    child = grid.boxes[pop!(parent.children)]
    Ljk = lcoeftrans(parent, child, p)
    # set the local coefs of the child due to colleague
    lcoef!(child,colparent,p,grid)
    println(Ljk)
    println(child.local_coef)

    num = 0
    den = 0
    for i = 1:length(Ljk)
        Ljkexpected = child.local_coef[i]
        num += abs(Ljkexpected - Ljk[i])^2
        den += abs(Ljkexpected)^2
    end

    E = sqrt(num/den)
    @test E < 10^-6

end

ltolTest1()
