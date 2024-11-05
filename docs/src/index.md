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
F(\tau) = F(-1) + \int_{-1}^{1}f(\tau)d\tau
```
where ``f(\tau)`` is an integrand that is valid on the domain [-1,1] and
``F(\tau)`` is an antiderivative of the integrand. 

!!! note
    If ``F(\tau)`` is an antiderivative of ``f(\tau)``, then ``F(\tau)+c`` is
    also an antiderivative of ``f(\tau)`` for any constant value ``c``. But, if
    ``F(-1)`` is known, then there is only one antiderivative that satisfies the
    integral.


## Methods

```@docs
interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)
clenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)
clenshaw_curtis_ivpi(N::Integer, M::Integer)
clenshaw_curtis_ivpii(N::Integer, M::Integer)
```