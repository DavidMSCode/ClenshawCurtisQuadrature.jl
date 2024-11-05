# ClenshawCurtisQuadrature.jl Documentation

## Summary

ClenshawCurtisQuadrature.jl is a Julia package that provides tools for
performing Clenshaw-Curtis quadrature, a numerical integration technique. This
method is particularly useful for quickly integrating arbitrary functions with
high accuracy. The package aims to offer efficient and easy-to-use functions
for integration tasks.

## Installation
ClenshawCurtisQuadrature.jl is currently in development and not included in any
package registries. Users will need to manually add the package with this
Github repository url
```Julia
julia> Using Pkg; Pkg.add(url="https://github.com/DavidMSCode/ClenshawCurtisQuadrature.jl")
```
or enter the package manager with 
```julia
julia> ]
```
and run
```julia
(@environment) pkg> add https://github.com/DavidMSCode/ClenshawCurtisQuadrature.jl
```

## Usage
The purpose of this module is to numerically solve integrals of the form

```math
F(\tau) = F(-1) + \int_{-1}^{\tau}f(\xi)d\xi
```
!!! note
    If ``F(\tau)`` is an antiderivative of ``f(\tau)``, then ``F(\tau)+c`` is
    also an antiderivative of ``f(\tau)`` for any constant value ``c``. But, if
    ``F(-1)`` is known, then there is only one antiderivative that satisfies the
    integral.

where ``f(\tau)`` is an integrand that is valid on the domain ``[-1,1]`` and
``F(\tau)`` is an antiderivative of the integrand. This is performed by
representing ``F(\tau)`` and ``f(\tau)`` as a summation of Chebyshev
polynomials on the domain ``[-1,1]``

```math
F(\tau) \approx \sum_{n=0}^N \alpha_n T_n(\tau)
```

```math
f(\tau) \approx \sum_{n=0}^{N-1} a_n T_n(\tau)
```
where ``T_n(\tau)`` is the ``n``-th Chebyshev polynomial of the first kind.
``N`` is the degree of the interpolating polynomial and ``a_n`` and ``\alpha_n``
are the interpolating coefficients of the integrand and the particular
antiderivative respectively. Because the integrals of the Chebyshev polynomials
are related the values of Chebyshev polynomials of higher and lower degrees, we
can perform linear operations on the integrand's coefficients to calculate
the corresponding antiderivative's coefficients. The functions in this module
are designed to generate the matrices required for this quadrature.

## Examples
If we want to know the integral on the domain ``[-1,1]`` of the function
``f(\tau)=\tau^2``, we can do so by first sampling the function at set of M+1
points. Because we know the antiderivative of ``f(\tau)`` is a 3rd order
polynomial, I'll choose ``N=3`` as our desired interpolating polynomial of the
solution. We can sample as many points as we want, but because the integrand
will be represented as a 2nd order polynomial we need at least 3 points to
calculate the coefficients.

To do so first generate the values of ``\tau`` at cosine spaced nodes.
!!! note
    Because we just wanted the value of the integral on the range ``[-1,1]`` 
    we assume that ``F(-1)=0``. The next examples will demonstrating solving
    an initial value problem when ``F(-1)\neq 0``.

    The Ta and T1 matrices returned by clenshaw\_curtis\_ivpi() are also
    interpolating matrices. They specifically contain the values of the
    chebyshev polynomials needed to interpolate the integrand and
    antiderivative at each of the ``M+1`` cosine spaced nodes we calculated
    before. Pre-calculating these matrices can be useful for iterative problems
    where the integrand is unkown as a function of ``\tau``.

```jldoctest
using ClenshawCurtisQuadrature

#set the degree variables
M = 2
N = 3

#Generate the cosine spaced nodes
Ms = 0:M
taus = cos.(pi*Ms/M)

#Sample the integrand
f = (x) -> x^2
ys = f.(taus)

#Get the Least Squares Operation and Quadrature matrix
A,Ta,P1,T1 = clenshaw_curtis_ivpi(N,M)

#Calculate the least squares fit interpolating polynomial coefficients
a = A*ys

#Apply the quadrature matrix
alpha = P1*a

#Get the interpolating polynomials at τ=1
T = interpolate(1,N)

F1 = T*alpha

# output

1-element Vector{Float64}:
 0.6666666666666666
```
The result is that ``\int_{-1}^1 \tau^2d\tau=\frac{2}{3}`` which is exactly the
correct answer.


## Methods

```@docs
interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)
clenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)
clenshaw_curtis_ivpi(N::Integer, M::Integer)
clenshaw_curtis_ivpii(N::Integer, M::Integer)
```