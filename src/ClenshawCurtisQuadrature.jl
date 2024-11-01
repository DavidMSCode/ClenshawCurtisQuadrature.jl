"""
Clenshaw Curtis Quadrature Module. Computes quadrature and interpolation matrices 
for performing Clenshaw-Curtis quadrature on the domain [-1,1].
"""
module ClenshawCurtisQuadrature

export clenshaw_curtis_nested_ivpd, interpolate, clenshaw_curtis_ivpi, clenshaw_curtis_ivpii

include("clenshawcurtisivp.jl")

end
