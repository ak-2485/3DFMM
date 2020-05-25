# An Adaptive FMM in Three Dimensions with the Basic Translation Operators

A Julia implementation of the adaptive 3D FMM described in *A Fast Adaptive Multipole Algorithm
in Three Dimensions* by Cheng, Greengard,and Rokhlin (1999), without the improvements to the translation and 
conversion operators. 

## Dependecies

In order to run a test, you will need the Julia packages GSL (for constructing the spherical harmonics) and LinearAlgebra. 

## Tests

You can run a test on the timing and accuracy of the implementation on a uniform random distribution of *num* particles in a domain *[-d,d]* with the maximum number of particles in a leaf *nmax* and expansion of order *p* using the following commands in Julia.

```
using FMMTests
FMMTests.uniformTest(num,nmax,p,d)
```
