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
    is_linear_stable(nw::Network, s0::NWState; marginally_stable=false, kwargs...)

Check if the fixpoint `s0` of the network `nw` is linearly stable by computing
the eigenvalues of the Jacobian matrix (or reduced Jacobian for constrained systems).

A fixpoint is linearly stable if all eigenvalues of the Jacobian have negative
real parts. For systems with algebraic constraints (non-identity mass matrix),
the reduced Jacobian is used following the approach in [1].
See [`jacobian_eigenvals`](@ref) for more details.

# Arguments
- `nw::Network`: The network dynamics object
- `s0::NWState`: The state to check for linear stability (must be a fixpoint)
- `marginally_stable::Bool=false`: If `true`, eigenvalues with zero real part are considered stable.
- `atol::Float64=1e-14`: Absolute tolerance for determining marginal stability. When `marginally_stable=true`, eigenvalues with `|real(λ)| < atol` are treated as having zero real part.
- `kwargs...`: Additional keyword arguments passed to `jacobian_eigenvals`

# Returns
- `Bool`: `true` if the fixpoint is linearly stable, `false` otherwise

# References
[1] "Power System Modelling and Scripting", F. Milano, Chapter 7.2.
"""
function is_linear_stable(nw::Network, s0; marginally_stable=false, atol=1e-14, kwargs...)
    isfixpoint(nw, s0) || error("The state s0 is not a fixpoint of the network nw.")
    λ = jacobian_eigenvals(nw, s0; kwargs...)
    comparator = if marginally_stable
        λ -> real(λ) - atol < 0.0 || isapprox(real(λ), 0.0; atol)
    else
        λ -> real(λ) < 0.0
    end
    if all(comparator, λ)
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
using LinearAlgebra

"""
    jacobian_participation(nw::Network, s0::NWState; normalize=true)

Compute participation factors showing which state variables participate in each mode.

# Arguments
- `nw::Network`: The network dynamics object
- `s0::NWState`: The state at which to compute participation factors
- `normalize=true`: Whether to normalize participation factors per mode to sum to 1

# Returns
- `eigvals`: Eigenvalues of the (reduced) Jacobian
- `pmat`: Participation matrix where pmat[k,i] is the participation of state k in mode i
- `d_idx`: Indices of differential states (relevant for constrained systems)

# Interpretation
For each mode i (column), pmat[:,i] shows which states are most involved.
Higher values indicate stronger participation.
"""
function jacobian_participation(nw::Network, s0; normalize=true)
    x0, p, t = uflat(s0), pflat(s0), s0.t
    M = nw.mass_matrix

    h!(dx, x) = nw(dx, x, p, t)
    j(x) = (dx = similar(x); ForwardDiff.jacobian(h!, dx, x))
    J = j(x0)

    if M == LinearAlgebra.I
        # Unconstrained system
        A_s = J
        d_idx = 1:size(J, 1)
    else
        # Constrained system - compute reduced Jacobian
        M != LinearAlgebra.Diagonal(M) && error("The constraints are not diagonal.")
        c_idx = findall(LinearAlgebra.diag(M) .== 0)
        d_idx = findall(LinearAlgebra.diag(M) .== 1)

        f_x = J[d_idx, d_idx]
        f_y = J[d_idx, c_idx]
        g_x = J[c_idx, d_idx]
        g_y = J[c_idx, c_idx]

        D = f_y * LinearAlgebra.pinv(g_y) * g_x
        A_s = f_x - D
    end

    # Compute eigenvalues and eigenvectors
    eigvals, V = eigen(A_s)  # Right eigenvectors (columns of V)
    W = inv(V)                # Left eigenvectors (rows of W)

    # Compute participation factors
    n_states = size(A_s, 1)
    n_modes = length(eigvals)
    pmat = zeros(n_states, n_modes)

    for i in 1:n_modes
        for k in 1:n_states
            pmat[k, i] = abs(W[i, k] * V[k, i])
        end
        if normalize
            pmat[:, i] ./= sum(pmat[:, i])  # Normalize so each mode sums to 1
        end
    end

    return eigvals, pmat, d_idx
end


"""
    show_mode_participation(nw::Network, s0::NWState, mode_idx::Int; threshold=0.05)

Display the states that significantly participate in a specific mode.

# Arguments
- `mode_idx`: Which mode/eigenvalue to analyze
- `threshold`: Only show states with participation > threshold (default: 0.05 = 5%)
"""
function show_mode_participation(nw::Network, s0, mode_idx::Int;
                                threshold=0.05)
    eigvals, pmat, d_idx = jacobian_participation(nw, s0)

    λ = eigvals[mode_idx]
    participations = pmat[:, mode_idx]

    println("Mode $mode_idx: λ = $λ")
    println("Frequency: $(abs(imag(λ))/(2π)) Hz, Damping: $(real(λ))")
    println("\nSignificant participants (> $(threshold*100)%):")

    state_names = string.(SII.variable_symbols(nw))
    for (i, p) in enumerate(participations)
        if p > threshold
            state_name = isnothing(state_names) ? "State $(d_idx[i])" : state_names[i]
            println("  $state_name: $(round(p*100, digits=2))%")
        end
    end
end


function cl_transfer_function(nw::Network, s0, input, output)
    # We assume the following system
    # M * ̇x = f(x, p)
    #     u = h(x, p) # input observable
    #     y = g(x, p) # output observable
    # we then linearize according to
    # M*δẋ = fx*δx + fp*δp
    #   δu = hx*δx + hp*δp
    #   δy = gx*δx + gp*δp
    # which is the following in freq domain
    # (M*s - fx)*x̂ = fp*p̂                          |x̂| = | Ms-fx  -fp |^(-1) |0| * û
    #            û = hx*x̂ + hp*p̂   together gives  |p̂|   |  hx     hp |      |1|
    #            ŷ = gx*x̂ + gp*p̂   which gives      ŷ = [gx px] * [x̂ p̂]ᵗ
    # therefore we get
    #  F(x) = [gx gu] * | Ms-fx  -fp |^(-1) |0|
    #                   |  hx     hp |      |1|

    M = nw.mass_matrix

    getu = SII.getu(nw, input)
    gety = SII.getu(nw, output)

    # first we need to define the jacobain hx, hp, gx, gp, fx, fp
    h = function(x, p)
        _s = NWState(nw, x, p, s0.t)
        getu(_s)
    end
    g = function(x, p)
        _s = NWState(nw, x, p, s0.t)
        gety(_s)
    end
    f = function(x, p)
        du = zeros(promote_type(eltype(x), eltype(p)), length(x))
        nw(du, x, p, s0.t)
        du
    end
    hx = ForwardDiff.jacobian(x -> h(x, pflat(s0)), uflat(s0))
    hp = ForwardDiff.jacobian(p -> h(uflat(s0), p), pflat(s0))
    gx = ForwardDiff.jacobian(x -> g(x, pflat(s0)), uflat(s0))
    gp = ForwardDiff.jacobian(p -> g(uflat(s0), p), pflat(s0))
    fx = ForwardDiff.jacobian(x -> f(x, pflat(s0)), uflat(s0))
    fp = ForwardDiff.jacobian(p -> f(uflat(s0), p), pflat(s0))

    dim_x = dim(nw)
    dim_u = length(input)

    F = function(s)
        X = [(M*s-fx)  -fp
            hx       hp]
        mask = [zeros(dim_x, dim_u)
                diagm(ones(dim_u))]
        [gx gp] * LinearAlgebra.pinv(X) * mask
    end
end

function cl_bus_tf_function(nw::Network, s0, sidx)
    i = resolvecompidx(nw, sidx)
    inputs = collect(VIndex(i, insym(nw[sidx])))
    # outputs = collect(VIndex(i, outsym(nw[sidx])))

    f = function(x, pert)
        du = zeros(promote_type(eltype(x), eltype(pert)), length(x))
        nw(du, x, pflat(s0), s0.t; input_perturb=Dict(VIndex(i) => pert))
        du
    end
    g = function(x, pert)
        du = f(x, pert) # make sure to fill buffers
        # XXX: Hacky way to get the cache
        outbuf = get_tmp(nw.caches.output, du) # get the outbuf which was filled by exec of the f function
        outbuf[nw.im.v_out[i]] # return the output of the bus component
    end

    M = nw.mass_matrix
    A = fx = ForwardDiff.jacobian(x -> f(x, zeros(length(inputs))), uflat(s0))
    B = fu = ForwardDiff.jacobian(pert -> f(uflat(s0), pert), zeros(length(inputs)))
    C = gx = ForwardDiff.jacobian(x -> g(x, zeros(length(inputs))), uflat(s0))
    D = gu = ForwardDiff.jacobian(pert -> g(uflat(s0), pert), zeros(length(inputs)))

    G(s) = C * ((M*s - A) \ B) + D
end

function cl_bus_invtf_function(nw::Network, s0, sidx)
    i = resolvecompidx(nw, sidx)
    # inputs = collect(VIndex(i, insym(nw[sidx])))
    outputs = collect(VIndex(i, outsym(nw[sidx])))

    f = function(x, pert)
        du = zeros(promote_type(eltype(x), eltype(pert)), length(x))
        nw(du, x, pflat(s0), s0.t; output_perturb=Dict(VIndex(i) => pert))
        du
    end
    g = function(x, pert)
        du = f(x, pert) # make sure to fill buffers
        # XXX: Hacky way to get the cache
        aggbuf = get_tmp(nw.caches.aggregation, du) # get the filled aggregation buffer (input for that vertex)
        aggbuf[nw.im.v_aggr[i]] # return the input of the bus component
    end

    M = nw.mass_matrix
    A = fx = ForwardDiff.jacobian(x -> f(x, zeros(length(outputs))), uflat(s0))
    B = fu = ForwardDiff.jacobian(pert -> f(uflat(s0), pert), zeros(length(outputs)))
    C = gx = ForwardDiff.jacobian(x -> g(x, zeros(length(outputs))), uflat(s0))
    D = gu = ForwardDiff.jacobian(pert -> g(uflat(s0), pert), zeros(length(outputs)))

    G(s) = C * ((M*s - A) \ B) + D
end

function cl_transfer_function2(nw::Network, s0; input, output)
    # output = collect(VIndex(i, insym(nw[sidx])))
    # input = collect(VIndex(i, outsym(nw[sidx])))
    M = nw.mass_matrix
    x0 = uflat(s0)
    p0 = pflat(s0)
    t0 = s0.t
    u0 = zeros(length(input))

    getu = SII.getu(nw, input)
    gety = SII.getu(nw, output)

    h = function(x)
        _s = NWState(nw, x, p0, t0)
        getu(_s)
    end
    hx = ForwardDiff.jacobian(h, uflat(s0))
    hx⁻¹ = LinearAlgebra.pinv(hx)

    f = function(x, u)
        x̂ = x + hx⁻¹*u
        du = zeros(eltype(x̂), length(x))
        nw(du, x̂, p0, t0)
        du
    end
    g = function(x, u)
        x̂ = x + hx⁻¹*u
        _s = NWState(nw, x̂, p0, t0)
        gety(_s)
    end

    A = fx = ForwardDiff.jacobian(x -> f(x, u0), x0)
    B = fu = ForwardDiff.jacobian(u -> f(x0, u), u0)
    C = gx = ForwardDiff.jacobian(x -> g(x, u0), x0)
    D = gu = ForwardDiff.jacobian(u -> g(x0, u), u0)

    G(s) = C * ((M*s - A) \ B) + D
end

function linearized_model(nw, s0; in, out)
    perturbation_channels = collect(in)
    observed_channels = collect(out)
    # perturbations need can only enter ant inputs or outputs of components (causal links)
    # Vertex Output -> pre gather, applyied to outbut
    # Vertex Input -> aggregation perturbation (needs to happen directly after loopback)
    # Edge Output -> applyed pre aggregate to outbuf
    # Edge Input -> requires buffered gathering, applyed post gather to gbuf
    δ0, perturb_maps = _classify_perturbation_channels(nw, perturbation_channels)
    M = nw.mass_matrix
    x0 = uflat(s0)
    p0 = pflat(s0)
    t0 = s0.t

    f = function(x, δ)
        du = zeros(promote_type(eltype(x), eltype(δ)), length(x))
        nw(du, x, p0, t0; perturb=δ, perturb_maps)
        du
    end
    obsf = SII.observed(nw, observed_channels)
    g = function(x, δ)
        out = zeros(promote_type(eltype(x), eltype(δ)), length(observed_channels))
        obsf(x, p0, t0, out; perturb=δ, perturb_maps)
        out
    end

    A = fx = ForwardDiff.jacobian(x -> f(x, δ0), x0)
    B = fu = ForwardDiff.jacobian(δ -> f(x0, δ), δ0)
    C = gx = ForwardDiff.jacobian(x -> g(x, δ0), x0)
    D = gu = ForwardDiff.jacobian(δ -> g(x0, δ), δ0)

    G(s) = C * ((M*s - A) \ B) + D
end

function _classify_perturbation_channels(nw, perturbation_channels)
    # we build a vector same size as perturbation
    δ0 = zeros(length(perturbation_channels))
    vo_map = Vector{NTuple{2, Int}}() # map outbuf idx -> δ idx
    vi_map = Vector{NTuple{2, Int}}() # map aggbuf idx -> δ idx
    eo_map = Vector{NTuple{2, Int}}() # map outbuf idx -> δ idx
    ei_map = Vector{NTuple{2, Int}}() # map gbuf idx -> δ idx

    # classify the perturbation channels
    for (i, sym) in enumerate(perturbation_channels)
        comp = getcomp(nw, sym)
        if comp isa VertexModel
            if sym.subidx ∈ outsym(comp)
                _idx = findfirst(==(sym.subidx), outsym(comp))
                outidx = nw.im.v_out[resolvecompidx(nw, sym)][_idx]
                push!(vo_map, (outidx, i))
            elseif sym.subidx ∈ insym(comp)
                _idx = findfirst(==(sym.subidx), insym(comp))
                aggidx = nw.im.v_aggr[resolvecompidx(nw, sym)][_idx]
                push!(vi_map, (aggidx, i))
            else
                throw(ArgumentError("Perturbation channel $sym is not an input or output of its component."))
            end
        elseif comp isa EdgeModel
            if sym.subidx ∈ outsym(comp).src
                _idx = findfirst(==(sym.subidx), outsym(comp).src)
                outidx = nw.im.e_out[resolvecompidx(nw, sym)].src[_idx]
                push!(eo_map, (outidx, i))
            elseif sym.subidx ∈ outsym(comp).dst
                _idx = findfirst(==(sym.subidx), outsym(comp).dst)
                outidx = nw.im.e_out[resolvecompidx(nw, sym)].dst[_idx]
                push!(eo_map, (outidx, i))
            elseif sym.subidx ∈ insym(comp).src
                _idx = findfirst(==(sym.subidx), insym(comp).src)
                gbufidx = nw.im.e_gbufr[resolvecompidx(nw, sym)].src[_idx]
                push!(ei_map, (gbufidx, i))
            elseif sym.subidx ∈ insym(comp).dst
                _idx = findfirst(==(sym.subidx), insym(comp).dst)
                gbufidx = nw.im.e_gbufr[resolvecompidx(nw, sym)].dst[_idx]
                push!(ei_map, (gbufidx, i))
            else
                throw(ArgumentError("Perturbation channel $sym is not an input or output of its component."))
            end
        else
            error()
        end
    end
    perturb_maps = (; vo_map, vi_map, eo_map, ei_map)
    δ0, perturb_maps
end
