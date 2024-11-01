using LinearAlgebra
"""
	interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)

Compute a matrix of Chebyshev polynomials of the first kind ``T_n(\\tau)`` for ``\\tau=taus`` and ``n = 0,1,...,N``.

# Arguments
- `taus::AbstractVector{<:Real}`: The points at which to evaluate the Chebyshev polynomials.
- `N::Integer`: The polynomial degree.
- `recursive::Bool`: If true, use the recursive formula to compute the Chebyshev polynomials within the domain [-1,1]. If false, use the trigonometric formulation. Default is false.

# Returns
- `Ts`: A matrix of Chebyshev polynomial values at the given values of tau.

# Description
This function computes the value of the N+1 Chebyshev polynomials at the given points in the domain [-1,1]. If the input `taus` is a scalar, it is converted to a vector. The function then computes the unweighted Chebyshev polynomial values at each tau.
"""
function interpolate(taus::AbstractVector{<:Real}, N::Integer; recursive::Bool = false)
	#find the indeces of the taus that are outside the domain
	idx_outside = findall(x -> abs(x) > 1, taus)
	if !isempty(idx_outside) #warn that if taus are outside domain then the output is extrapolation and has no guarantee of accuracy
		@warn ("The input values at the following indices are outside the domain [-1,1]: $idx_outside. Extrapolated values may be innacurate.")
	end

	#Get the uneweighted Chebyshev polynomial values at each tau
	js = 0:N
	if !recursive
		#Default behaviour: Use the trig formulation and only use the recursive formula if tau is outside the domain [-1,1]
		Ts = [abs(tau) <= 1 ? trig_chebyshev(tau, j) : recursive_chebyshev(tau, j) for tau in taus, j in js]
	else
		#always use the recursive formula if recursive is true
		Ts = [recursive_chebyshev(tau, j) for tau in taus, j in js]
	end
	return Ts
end

"""
	interpolate(tau::Real, N::Integer; recursive = false)

Computes a 1x(N+1) matrix of Chebyshev polynomials of the first kind ``T_n(\\tau)`` for ``n = 0,1,...,N``.

# Arguments
- `tau::Real`: The point at which to evaluate the Chebyshev polynomials.
- `N::Integer`: The polynomial degree.
- `recursive::Bool`: If true, use the recursive formula to compute the Chebyshev polynomials within the domain [-1,1]. If false, use the trigonometric formulation. Default is false.

# Returns
- `Ts`: A matrix of Chebyshev polynomial values at the given value of tau.

# Description
This function computes the value of the N+1 Chebyshev polynomials at the given point in the domain [-1,1].  The function then computes the unweighted Chebyshev polynomial at tau.
"""
function interpolate(tau::Real, N::Integer; recursive = false)
	return interpolate([tau], N; recursive = recursive)
end

"""
	trig_chebyshev(tau::Real, N::Integer)

Computes cos(N * acos(tau)). The trigonometric form of the Chebyshev polynomial N at the point tau.

# Arguments
- `tau::Real`: The point at which to evaluate the Chebyshev polynomial.
- `N::Integer`: The polynomial degree.

# Returns
- `T`: The Chebyshev polynomial value at the given point.
"""
function trig_chebyshev(tau::Real, N::Integer)
	return cos(N * acos(tau))
end

"""
	recursive_chebyshev(tau::Real, N::Integer)

Computes ``T_N(tau) = 2 * tau * T_{N-1}(tau) - T_{N-2}(tau)`` the Nth Chebyshev polynomial at tau using the recursive formula.

# Arguments
- `tau::Real`: The point at which to evaluate the Chebyshev polynomial.
- `N::Integer`: The polynomial degree.

# Returns
- `T`: The Chebyshev polynomial value at the given point.
"""
function recursive_chebyshev(tau::Real, N::Integer)
	if N == 0
		return 1
	elseif N == 1
		return tau
	else
		return 2 * tau * recursive_chebyshev(tau, N - 1) - recursive_chebyshev(tau, N - 2)
	end
end

"""
	chebyshev(N::Integer, M::Integer)

Compute the value of the N+1 Chebyshev polynomials at the M+1 points.

# Arguments
- `N::Integer`: The polynomial degree.
- `M::Integer`: The sampling degree. Must be greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.

# Returns
- `Ts`: The Chebyshev polynomial values at the M+1 points.

# Description
This function computes the value of the N+1 Chebyshev polynomials at the M+1 cosine spaced nodes. It first generates the cosine spaced sample points from [-1,1] using the `cosineSamplePoints` function. Then, it computes the unweighted Chebyshev polynomial values at each tau.
"""
function chebyshev(N::Integer, M::Integer)
	# Get the uneweighted Chebyshev polynomials at each of the M+1 points
	taus = cosineSamplePoints(M)
	js = 0:N
	Ts = interpolate(taus, N)
	return Ts
end

"""
	cosineSamplePoints(M::Integer)

Compute the cosine spaced sample points from [-1,1].

# Arguments
- `M::Integer`: The number of sample points.

# Returns
- `taus`: The cosine spaced sample points from [-1,1].

# Description
This function computes the cosine spaced sample points from [-1,1]. It generates M+1 points from 0 to M and computes the cosine of each point multiplied by pi divided by M.
"""
function cosineSamplePoints(M::Integer)
	# Cosine Sample points (M+1 points) from [-1,1]
	taus = 0:M
	return -cos.(pi * taus / M)
end

"""
	lsq_chebyshev_fit(N::Integer, M::Integer)

Compute the Chebyshev polynomials and the Least Squares Operator matrix.

# Arguments
- `N::Integer`: The polynomial degree.
- `M::Integer`: The sampling degree. Must be greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.

# Returns
- `T`: The Chebyshev polynomial values at the M+1 points.
- `A`: The Least Squares Operator matrix.

# Description
This function computes the Chebyshev polynomials at the M+1 points and the Least Squares Operator matrix. It first generates the Chebyshev polynomials using the `chebyshev` function. Then, it constructs the weights matrix, W, and the V matrix. Finally, it computes the Least Squares Operator matrix, A.

# Example
```julia
N = 5
M = 5
T, A = lsq_chebyshev_fit(N, M)
```
"""
function lsq_chebyshev_fit(N::Integer, M::Integer)
	#generate the Chebyshev polynomials at the nodes
	T = chebyshev(N, M)

	#weights matrix
	wdiag = [0.5; ones(M - 1); 0.5]
	W = Diagonal(wdiag)

	#V matrix
	vdiag = [1 / M; 2 .* ones(N) ./ M]
	if M == N
		vdiag[end] = 1 / M
	end
	V = Diagonal(vdiag)

	#transpose T
	Tt = transpose(T)

	#Least Squares Operator
	A = V * Tt * W

	return T, A
end

"""
	clenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)

Compute the Clenshaw-Curtis quadrature and Chebyshev basis function matrices for a d-th order integral.

# Arguments
- `N::Integer`: The polynomial degree.
- `M::Integer`: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.
- `d::Integer`: The integral order.

# Returns
- `A`: The Least Squares Operator matrix.
- `P`: The Quadrature Matrix.
- `T`: The Chebyshev Matrix.

# Description
This function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the `lsq_chebyshev_fit` function. Then, it calculates the Constants of Integration, and constructs the S matrix. Finally, it computes the Clenshaw Curtis Quadrature matrix.
"""
function clenshaw_curtis_nested_ivpd(N::Integer, M::Integer, d::Integer)
	if M < N
		throw(ArgumentError("The number of sampling nodes must be greater than the polynomial order, N."))
	end
	if d > N
		throw(ArgumentError("The polynomial order N must be greater than or equal to the integral order d."))
	end

	# Least Squares Operator for "acceleration"
	Ta, A = lsq_chebyshev_fit(N - d, M)

	# Constants of Integration
	ks = 0:N
	Lrow = cos.(ks * pi)
	L = vcat(Lrow', zeros(N, N + 1))

	#S matrix
	tempdiag = [1; [1 / (2 * i) for i in 1:N]]
	temp = Diagonal(tempdiag)
	temp2 = diagm(N + 1, N, -1 => ones(N), 1 => [0; -ones(N - 2)])
	S = temp * temp2
	S[1, 1] = 0.25
	S[2, 1] = 1.0

	#Integration Operator
	temp3 = -L + Diagonal(ones(N + 1))
	P = temp3 * S

	# Chebyshev Matrix
	T = chebyshev(N, M)

	return A, P, T
	return
end

"""
	clenshaw_curtis_ivpii(N::Integer, M::Integer)

Compute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a second order initial value problem.

# Arguments
- `N::Integer`: The polynomial degree.
- `M::Integer`: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.

# Returns
- `A`: The Least Squares Operator matrix.
- `Ta`: The "acceleration" Chebyshev Matrix.
- `P1`: The Quadrature Matrix for acceleration to velocity.
- `T1`: The "Velocity" Chebyshev Matrix.
- `P2`: The Quadrature Matrix for velocity to position.
- `T2`: The "Position" Chebyshev Matrix.

# Description
This function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the `lsq_chebyshev_fit` function. Then, it calculates the "Position" and "Velocity" Constants of Integration, and constructs the S matrices for "velocity" and "position". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity and velocity to position.

# Example
```julia
N = 5
M = 5
A, Ta, P1, T1, P2, T2 = clenshaw_curtis_ivpii(N, M)
"""
function clenshaw_curtis_ivpii(N::Integer, M::Integer)
	if M < N
		throw(ArgumentError("The number of sampling nodes must be greater than the polynomial order, N."))
	end

	# Least Squares Operator for "acceleration"
	Ta, A = lsq_chebyshev_fit(N - 2, M)

	# "Position" Constants of Integration
	ks = 0:N
	Lprow = cos.(ks * pi)
	Lp = vcat(Lprow', zeros(N, N + 1))

	#S matrix for "position"
	temp4diag = [1; [1 / (2 * i) for i in 1:N]]
	temp4 = Diagonal(temp4diag)
	temp5 = diagm(N + 1, N, -1 => ones(N), 1 => [0; -ones(N - 2)])
	Sp = temp4 * temp5
	Sp[1, 1] = 0.25
	Sp[2, 1] = 1.0

	#Picard Integration Operator for velocity to position
	temp6 = -Lp + Diagonal(ones(N + 1))
	P2 = temp6 * Sp

	# "Velocity" Quadrature matrix is subset of the "Position" Quadrature matrix
	P1 = P2[1:end-1, 1:end-1]

	# "Position" Chebyshev Matrix
	T2 = chebyshev(N, M)
	# "Velocity" Chebyshev Matrix is subset of the "Position" Chebyshev Matrix
	T1 = T2[:, 1:end-1]

	return A, Ta, P1, T1, P2, T2
end

"""
	clenshaw_curtis_ivpi(N::Integer, M::Integer)

Compute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a first order initial value problem.

# Arguments
- `N::Integer`: The polynomial degree.
- `M::Integer`: The sampling degree. Must bes greater than or equal to the polynomial degree. This is equal to the total number of function sampling points minus 1.

# Returns
- `A`: The Least Squares Operator matrix.
- `Ta`: The "acceleration" Chebyshev Matrix.
- `P1`: The Quadrature Matrix for acceleration to velocity.
- `T1`: The "Velocity" Chebyshev Matrix.

# Description
This function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the `lsq_chebyshev_fit` function. Then, it calculates the "Velocity" Constants of Integration, and constructs the S matrices for "velocity". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity.
"""
function clenshaw_curtis_ivpi(N::Integer, M::Integer)
	if M < N
		throw(ArgumentError("The number of sampling nodes, M, must be greater than or equal to the polynomial order, N."))
	end

	# Least Squares Operator for "acceleration"
	Ta, A = lsq_chebyshev_fit(N - 1, M)

	# "Velocity" Constants of Integration
	ks = 0:N
	Lvrow = cos.(ks * pi)
	Lv = vcat(Lvrow', zeros(N, N + 1))

	#S matrix for "velocity"
	temp1diag = [1; [1 / (2 * i) for i in 1:N]]
	temp1 = Diagonal(temp1diag)
	temp2 = diagm(N + 1, N, -1 => ones(N), 1 => [0; -ones(N - 2)])
	Sv = temp1 * temp2
	Sv[1, 1] = 0.25
	Sv[2, 1] = 1.0

	#Picard Integration Operator for accleration to velocity
	temp3 = -Lv + Diagonal(ones(N + 1))
	P1 = temp3 * Sv

	# "Velocity" Chebyshev Matrix
	T1 = chebyshev(N, M)

	return A, Ta, P1, T1
end

"""
	clenshaw_curtis_ivpii(N::Integer)

Compute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a second order initial value problem.

# Arguments
- `N::Integer`: The polynomial degree.

# Returns
- `A`: The Least Squares Operator matrix.
- `Ta`: The "acceleration" Chebyshev Matrix.
- `P1`: The Quadrature Matrix for acceleration to velocity.
- `T1`: The "Velocity" Chebyshev Matrix.
- `P2`: The Quadrature Matrix for velocity to position.
- `T2`: The "Position" Chebyshev Matrix.

# Description
This function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the `lsq_chebyshev_fit` function. Then, it calculates the "Position" and "Velocity" Constants of Integration, and constructs the S matrices for "velocity" and "position". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity and velocity to position.

"""
function clenshaw_curtis_ivpii(N::Integer)
	return clenshaw_curtis_ivpii(N, N)
end

"""
	clenshaw_curtis_ivpi(N::Integer)

Compute the Clenshaw-Curtis quadrature and Cebyshev basis function matrices for a first order initial value problem.

# Arguments
- `N::Integer`: The polynomial degree.

# Returns
- `A`: The Least Squares Operator matrix.
- `Ta`: The "acceleration" Chebyshev Matrix.
- `P1`: The Quadrature Matrix for acceleration to velocity.
- `T1`: The "Velocity" Chebyshev Matrix.

# Description
This function computes the Clenshaw-Curtis quadrature matrices and the basis function vectors a. It first generates the Chebyshev polynomials and the Least Squares Operator matrix using the `lsq_chebyshev_fit` function. Then, it calculates the "Velocity" Constants of Integration, and constructs the S matrices for "velocity". Finally, it computes the Clenshaw Curtis Quadrature matrices for acceleration to velocity.
"""
function clenshaw_curtis_ivpi(N::Integer)
	return clenshaw_curtis_ivpi(N, N)
end