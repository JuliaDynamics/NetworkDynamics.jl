"""
    isfixpoint(nw::Network, s0::NWState; tol=1e-10)

Check if the state `s0` is a fixpoint of the network `nw` by
calculating the the RHS and check that every entry is within the
given tolerance `tol`.
"""
function isfixpoint(nw::Network, s0::NWState; tol=1e-10)
    # Check if the state is a fixpoint of the network
    u0 = uflat(s0)
    p0 = pflat(s0)
    du = zeros(eltype(u0), length(u0))
    nw(du, u0, p0, s0.t)

    # Check if the change is within the tolerance
    return all(abs.(du) .< tol)
end

"""
    is_linear_stable(nw::Network, s0::NWState; kwargs...)

Check if the fixpoint `s0` of the network `nw` is linearly stable by computing
the eigenvalues of the Jacobian matrix (or reduced Jacobian for constrained systems).

A fixpoint is linearly stable if all eigenvalues of the Jacobian have negative
real parts. For systems with algebraic constraints (non-identity mass matrix),
the reduced Jacobian is used following the approach in [1].
See [`jacobian_eigenvals`](@ref) for more details.

# Arguments
- `nw::Network`: The network dynamics object
- `s0::NWState`: The state to check for linear stability (must be a fixpoint)
- `kwargs...`: Additional keyword arguments passed to `jacobian_eigenvals`

# Returns
- `Bool`: `true` if the fixpoint is linearly stable, `false` otherwise

# References
[1] "Power System Modelling and Scripting", F. Milano, Chapter 7.2.
"""
function is_linear_stable(nw::Network, s0; kwargs...)
    isfixpoint(nw, s0) || error("The state s0 is not a fixpoint of the network nw.")
    λ = jacobian_eigenvals(nw, s0; kwargs...)
    if all(λ -> real(λ) < 0.0, λ)
        return true
    else
        return false
    end
end

"""
    jacobian_eigenvals(nw::Network, s0::NWState; eigvalf=LinearAlgebra.eigvals)

Compute the eigenvalues of the Jacobian matrix for linear stability analysis of
the network dynamics at state `s0`.

For systems without algebraic constraints (identity mass matrix), this returns
the eigenvalues of the full Jacobian matrix. For constrained systems (non-identity
mass matrix), it computes the eigenvalues of the reduced Jacobian following the
approach for differential-algebraic equations outlined in [1]

# Arguments
- `nw::Network`: The network dynamics object
- `s0::NWState`: The state at which to compute the Jacobian eigenvalues
- `eigvalf`: Function to compute eigenvalues (default: `LinearAlgebra.eigvals`)

# Returns
- `Vector`: Eigenvalues of the Jacobian (or reduced Jacobian for constrained systems)

# Algorithm
For unconstrained systems (M = I):
- Computes eigenvalues of the full Jacobian J

For constrained systems (M ≠ I, differential-algebraic equations):
- The system has the form: M * dz/dt = f(z, t) where M is a diagonal mass matrix
- Variables are partitioned into differential (M_ii = 1) and algebraic (M_ii = 0) components
- Let z = [x; y] where x are differential and y are algebraic variables
- The Jacobian J = ∂f/∂z is partitioned as:
  ```
  J = [f_x  f_y]  where f_x = ∂f_d/∂x, f_y = ∂f_d/∂y
      [g_x  g_y]        g_x = ∂g_a/∂x, g_y = ∂g_a/∂y
  ```
- For the algebraic constraints 0 = g_a(x, y), we have dy/dt = -g_y^(-1) * g_x * dx/dt
- Substituting into the differential equations gives the reduced system:
  dx/dt = (f_x - f_y * g_y^(-1) * g_x) * x = A_s * x
- The eigenvalues of the reduced Jacobian A_s determine stability
- This approach follows the theory of differential-algebraic equations [1]

# References
[1] "Power System Modelling and Scripting", F. Milano, Chapter 7.2.
"""
function jacobian_eigenvals(nw::Network, s0; eigvalf=LinearAlgebra.eigvals)
    x0, p, t = uflat(s0), pflat(s0), s0.t # unpack state
    M = nw.mass_matrix

    h!(dx, x) = nw(dx, x, p, t)
    j(x) = (dx = similar(x); ForwardDiff.jacobian(h!, dx, x)) # Jacobian
    J = j(x0) # Full Jacobian at the equilibrium point

    if M == LinearAlgebra.I # Constraint free system -> use eigenvalues of jacobian
        return  eigvalf(J)
    else # Constraints -> Use eigenvalues of reduced jacobian
        M != LinearAlgebra.Diagonal(M) && error("The constraints are not diagonal.")
        c_idx = findall(LinearAlgebra.diag(M) .== 0)
        d_idx = findall(LinearAlgebra.diag(M) .== 1)

        f_x = J[d_idx, d_idx] # Differential equations evaluated at the differential variables
        f_y = J[d_idx, c_idx] # Differential equations evaluated at the constrained variables

        g_x = J[c_idx, d_idx] # Constrained equations evaluated at the differential variables
        g_y = J[c_idx, c_idx] # Constrained equations evaluated at the constrained variables

        D = f_y * LinearAlgebra.pinv(g_y) * g_x # Degradation matrix
        A_s = f_x - D             # State matrix / Reduced Jacobian (eq. 7.16 in [1])
        return eigvalf(A_s)       # Eigenvalues of the reduced jacobian
    end
end
