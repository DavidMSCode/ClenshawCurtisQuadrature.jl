var documenterSearchIndex = {"docs":
[{"location":"#ClenshawCurtisQuadrature.jl-Documentation","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.jl Documentation","text":"","category":"section"},{"location":"","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.jl Documentation","text":"interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)\nclenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)\nclenshaw_curtis_ivpi(N::Integer, M::Integer)\nclenshaw_curtis_ivpii(N::Integer, M::Integer)","category":"page"},{"location":"#ClenshawCurtisQuadrature.interpolate-Tuple{AbstractVector{<:Real}, Integer}","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.interpolate","text":"interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)\n\nCompute a matrix of Chebyshev polynomials of the first kind T_n(tau) for tau=taus and n = 01N.\n\nArguments\n\ntaus::AbstractVector{<:Real}: The points at which to evaluate the Chebyshev polynomials.\nN::Integer: The polynomial degree.\nrecursive::Bool: If true, use the recursive formula to compute the Chebyshev polynomials within the domain [-1,1]. If false, use the trigonometric formulation. Default is false.\n\nReturns\n\nTs: A matrix of Chebyshev polynomial values at the given values of tau.\n\nDescription\n\nThis function computes the value of the N+1 Chebyshev polynomials at the given points in the domain [-1,1]. If the input taus is a scalar, it is converted to a vector. The function then computes the unweighted Chebyshev polynomial values at each tau.\n\n\n\n\n\n","category":"method"},{"location":"#ClenshawCurtisQuadrature.clenshaw_curtis_nested_ivpd-Tuple{Integer, Integer, Integer}","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.clenshaw_curtis_nested_ivpd","text":"clenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)\n\nCompute the Clenshaw-Curtis quadrature and Chebyshev basis function matrices for a d-th order integral.\n\nArguments\n\nN::Integer: The polynomial degree.\nM::Integer: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.\nd::Integer: The integral order.\n\nReturns\n\nA: The Least Squares Operator matrix.\nP: The Quadrature Matrix.\nT: The Chebyshev Matrix.\n\nDescription\n\nThis function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the lsq_chebyshev_fit function. Then, it calculates the Constants of Integration, and constructs the S matrix. Finally, it computes the Clenshaw Curtis Quadrature matrix.\n\n\n\n\n\n","category":"method"},{"location":"#ClenshawCurtisQuadrature.clenshaw_curtis_ivpi-Tuple{Integer, Integer}","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.clenshaw_curtis_ivpi","text":"clenshaw_curtis_ivpi(N::Integer, M::Integer)\n\nCompute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a first order initial value problem.\n\nArguments\n\nN::Integer: The polynomial degree.\nM::Integer: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.\n\nReturns\n\nA: The Least Squares Operator matrix.\nTa: The \"acceleration\" Chebyshev Matrix.\nP1: The Quadrature Matrix for acceleration to velocity.\nT1: The \"Velocity\" Chebyshev Matrix.\n\nDescription\n\nThis function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the lsq_chebyshev_fit function. Then, it calculates the \"Velocity\" Constants of Integration, and constructs the S matrices for \"velocity\". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity.\n\n\n\n\n\n","category":"method"},{"location":"#ClenshawCurtisQuadrature.clenshaw_curtis_ivpii-Tuple{Integer, Integer}","page":"ClenshawCurtisQuadrature.jl Documentation","title":"ClenshawCurtisQuadrature.clenshaw_curtis_ivpii","text":"clenshaw_curtis_ivpii(N::Integer, M::Integer)\n\nCompute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a second order initial value problem.\n\nArguments\n\nN::Integer: The polynomial degree.\nM::Integer: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.\n\nReturns\n\nA: The Least Squares Operator matrix.\nTa: The \"acceleration\" Chebyshev Matrix.\nP1: The Quadrature Matrix for acceleration to velocity.\nT1: The \"Velocity\" Chebyshev Matrix.\nP2: The Quadrature Matrix for velocity to position.\nT2: The \"Position\" Chebyshev Matrix.\n\nDescription\n\nThis function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the lsq_chebyshev_fit function. Then, it calculates the \"Position\" and \"Velocity\" Constants of Integration, and constructs the S matrices for \"velocity\" and \"position\". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity and velocity to position.\n\nExample\n\n```julia N = 5 M = 5 A, Ta, P1, T1, P2, T2 = clenshawcurtisivpii(N, M)\n\n\n\n\n\n","category":"method"}]
}