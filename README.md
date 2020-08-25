# An adaptive FMM in three dimensions with basic translation operators

A Julia implementation of the adaptive 3D FMM described in *A Fast Adaptive Multipole Algorithm
in Three Dimensions* by Cheng, Greengard,and Rokhlin (1999), without the improvements to the translation and 
conversion operators. 

## Dependencies

In order to run an example or test, you will need the Julia packages GSL (for constructing the spherical harmonics), LinearAlgebra, Random, and Distributions. 

## Examples
### Calculate the potential at each particle in a uniform random distribution of charges
In order to generate a uniform random distribution of of *num* particles in a domain *[-d,d]*, type the following at the Julia prompt.

```
using EnsembleTests
particles, _,minbound,maxbound = rndmparticledist3(num::Int,d::Float)
```

Then, to apply the FMM to the distribution with the maximum number of particles in a leaf *nmax* and expansion of order *p*, type the following at the Julia prompt.

```
using FMMTests
ϕ = fastmultipole(particles,minbound,maxbound,nmax::Int,p::Int)
```

## Tests

You can run a test on the timing and accuracy of the implementation on a uniform random distribution of *num* particles in a domain *[-d,d]* with the maximum number of particles in a leaf *nmax* and expansion of order *p* using the following commands in Julia.

```
import FMMTests
FMMTests.uniformTest(num::Int,nmax::Int,p::Int,d::Float)
```
