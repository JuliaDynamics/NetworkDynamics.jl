"""
    NetworkDescriptorSystem

Descriptor system representation

    M ẋ = A x + B u
      y = C x + D u

of a linearized network. Returned by [`linearize_network`](@ref).

The struct is callable as a transfer function: `sys(s)` computes
`C (M s - A)⁻¹ B + D`. For SISO systems (scalar `in_syms`/`out_syms`) the
result is a scalar.

Access the properties using `sys.A`, `.B`, `.C`, `.D`, `.M` to obtain matrices.
Use `sym(⋅)`, `insym(⋅)` and `outsym(⋅)` to obtain state names (`sym` is
optional and may be `nothing`).

See also: [`linearize_network`](@ref), [`reduce_dae`](@ref)
"""
@kwdef struct NetworkDescriptorSystem{MT,AT,BT,CT,DT,IST,OST}
    M::MT = LinearAlgebra.I
    A::AT
    B::BT = nothing
    C::CT = nothing
    D::DT = nothing
    sym::Union{Vector, Nothing} = nothing
    insym::IST = nothing
    outsym::OST = nothing
    function NetworkDescriptorSystem(M, A, B, C, D, sym, insym, outsym)
        dim = size(A, 1)
        if !isnothing(sym) && length(sym) != dim
            throw(ArgumentError("Length of sym must match number of rows in A"))
        end
        if any(isnothing, (B, C, D, insym, outsym)) && !all(isnothing, (B, C, D, insym, outsym))
            throw(ArgumentError("B, C, D, insym, and outsym must all be provided or all be nothing"))
        end
        if !isnothing(B) # has in/out
            if size(B, 1) != dim
                throw(ArgumentError("Number of rows in B must match number of rows in A"))
            end
            if size(C, 2) != dim
                throw(ArgumentError("Number of columns in C must match number of rows in A"))
            end
            if size(D, 1) != size(C, 1) || size(D, 2) != size(B, 2)
                throw(ArgumentError("D must have same number of rows as C and same number of columns as B"))
            end
            _nsym(s::SymbolicIndex) = 1
            _nsym(s) = length(s)
            if _nsym(insym) != size(B, 2)
                throw(ArgumentError("Length of insym must match number of columns in B"))
            end
            if _nsym(outsym) != size(C, 1)
                throw(ArgumentError("Length of outsym must match number of rows in C"))
            end
        end
        if M != LinearAlgebra.I
            if M isa UniformScaling
                if M.λ != 1
                    throw(ArgumentError("M must be identity or a uniform scaling of 1"))
                end
            else # matrix
                if !(LinearAlgebra.isdiag(M))
                    throw(ArgumentError("M must be diagonal!"))
                end
                if any(x -> x != 0 && x != 1, LinearAlgebra.diag(M))
                    throw(ArgumentError("M must be a mass matrix with diagonal entries of 0 or 1"))
                end
                if size(M, 1) != dim || size(M, 2) != dim
                    throw(ArgumentError("M must be square and match number of rows in A"))
                end
            end
        end
        new{typeof(M), typeof(A), typeof(B), typeof(C), typeof(D), typeof(insym), typeof(outsym)}(M, A, B, C, D, sym, insym, outsym)
    end
end

function (sys::NetworkDescriptorSystem)(s)
    val = sys.C * ((sys.M * s - sys.A) \ sys.B) + sys.D
    if sys.insym isa SymbolicIndex && sys.outsym isa SymbolicIndex
        return only(val)
    end
    return val
end

dim(sys::NetworkDescriptorSystem) = size(sys.A, 1)
indim(sys::NetworkDescriptorSystem) = sys.B === nothing ? 0 : size(sys.B, 2)
outdim(sys::NetworkDescriptorSystem) = sys.C === nothing ? 0 : size(sys.C, 1)
sym(sys::NetworkDescriptorSystem) = sys.sym
insym(sys::NetworkDescriptorSystem) = sys.insym
outsym(sys::NetworkDescriptorSystem) = sys.outsym
constraint_dim(sys::NetworkDescriptorSystem{<:UniformScaling}) = 0
constraint_dim(sys::NetworkDescriptorSystem) = dim(sys) - sum(sys.M)

function Base.show(io::IO, ::MIME"text/plain", sys::NetworkDescriptorSystem)
    compact = get(io, :compact, false)::Bool
    printstyled(io, "NetworkDescriptorSystem"; bold=true)
    if sys.B === nothing
        if compact
            print(io, "(dim=$(dim(sys)))")
        else
            print(io, " with\n")
            print(io, "  ├─ $(dim(sys)) states")
            nc = constraint_dim(sys)
            @assert nc >= 0
            if nc == 0
                print(io, " (no algebraic constraints)\n")
            else
                constraints = maybe_plural(dim(sys) - sum(sys.M), "constaint")[2]
                print(io, " ($(nc) algebraic $(constraints))\n")
            end
            print(io, "  └─ without inputs/outputs")
        end
    else
        if compact
            print(io, "(dim=$(dim(sys)), in=$(indim(sys)), out=$(outdim(sys)))")
        else
            print(io, " with")
            strvec = String[]
            push!(strvec, "$(dim(sys)) &states")
            ni, _inputs = maybe_plural(indim(sys), "input")
            no, _outputs = maybe_plural(outdim(sys), "output")
            push!(strvec, "$(ni) &$(_inputs): &$(insym(sys))")
            push!(strvec, "$(no) &$(_outputs): &$(outsym(sys))")
            print_treelike(io, align_strings(strvec))
        end
    end
end

"""
    isfixpoint(s0::NWState; tol=1e-10)

Check if the state `s0` is a fixpoint of the network by
calculating the the RHS and check that every entry is within the
given tolerance `tol`.
"""
function isfixpoint(s0::NWState; tol=1e-10)
    nw = extract_nw(s0)
    u0 = uflat(s0)
    p0 = pflat(s0)
    du = zeros(eltype(u0), length(u0))
    nw(du, u0, p0, s0.t)
    return all(abs.(du) .< tol)
end

"""
    linearize_network(s0::NWState)
    linearize_network(s0::NWState; in, out)

Linearize the network dynamics around state `s0`.

Without `in`/`out`, returns a [`NetworkDescriptorSystem`](@ref) with `M`, `A`,
and `sym`.

With `in`/`out`, returns a `NetworkDescriptorSystem` with full state-space
matrices (`M`, `A`, `B`, `C`, `D`) for the perturbation channels specified by
`in` and observations specified by `out`. The returned system is callable as a
transfer function: `sys(s) = C * (M*s - A)⁻¹ * B + D`.

The `in` keyword agument specifies the perturbance channel.
Perturbance channels are related to the input/output structure of the components
(input indices musst match parameters or input or output states of individual components).
They can be either
- parameters (where `δp` is additive perturbation)
- `δφ_in (edge)`: Perturbation of the edge input (enters just this edge)
- `δφ_out (vertex)`: Perturbation of the vertex output (enters all connected edges)
- `δΦ_out (edge)` and `δΦ_in (vertex)`: Perturbation of the flow before or after aggregation. Mathematicially equivalent, only enters the vertex anyway.

```
                  δφ_in          more edges         δφ_in
                    ↓                △                ↓
                    ∙────────────────┼────────────────∙
n ⋯───╮             │                │ potential      │             ╭───⋯ n
e     │             │        δφ_out →∙  φ out         │             │     e
x  ┏━━▽━━━━━━━━━━━━━▽━━┓   ╔═════════△═════════╗   ┏━━▽━━━━━━━━━━━━━▽━━┓  x
t  ┃ EdgeModel         ┃   ║ VertexModel       ║   ┃ EdgeModel         ┃  t
   ┃ ẋ = f(x, φ, p, t) ┃   ║ ẋ = f(x, Φ, p, t) ║   ┃ ẋ = f(x, φ, p, t) ┃
n  ┃ Φ = g(x, φ, p, t) ┃   ║ φ = g(x, p, t)    ║   ┃ Φ = g(x, φ, p, t) ┃  n
o  ┗━━▽━━━━━━━━━━━━━▽━━┛   ╚═════════△═════════╝   ┗━━▽━━━━━━━━━━━━━▽━━┛  o
d     │        flow │ Φ out   δΦ_in →∙           flow │ Φ out       │     d
e ⋯───╯             │              ╭─┴─╮              │             ╰───⋯ e
            δΦ_out →∙──────────────▷ + ◁──────────────∙← δΦ_out
                                   ╰─△─╯
                                     │
                                 more edges
```

The output channel `out` can by any observable of the system.

# Examples
```julia
# Simple linearization for eigenvalue analysis
sys = linearize_network(s0)
eigvals = jacobian_eigenvals(sys)

# Full ABCD linearization with transfer function
sys = linearize_network(s0; in=[VIndex(1, :u_r)], out=[VIndex(1, :i_r)])
sys(im * 2π * 50)  # evaluate transfer function at 50 Hz
```

See also: [`jacobian_eigenvals`](@ref), [`participation_factors`](@ref)
"""
function linearize_network(s0::NWState; in=nothing, out=nothing)
    nw = extract_nw(s0)
    x0 = uflat(s0)
    p0 = pflat(s0)
    t0 = s0.t
    M = nw.mass_matrix
    sym = SII.variable_symbols(nw)

    if isnothing(in) && isnothing(out)
        h!(dx, x) = nw(dx, x, p0, t0)
        A = ForwardDiff.jacobian(h!, similar(x0), x0)
        return NetworkDescriptorSystem(; M, A, sym)
    else
        isnothing(in) && throw(ArgumentError("Must specify `in` when `out` is given"))
        isnothing(out) && throw(ArgumentError("Must specify `out` when `in` is given"))

        # check for siso
        single_in = in isa SymbolicIndex && SII.symbolic_type(in) == SII.ScalarSymbolic()
        single_out = out isa SymbolicIndex && SII.symbolic_type(out) == SII.ScalarSymbolic()
        siso = single_in && single_out

        insym = single_in ? [in] : collect(in)
        outsym = single_out ? [out] : collect(out)
        δ0, perturb_maps = _classify_perturbation_channels(nw, insym)

        _pdim = pdim(nw)
        f = function(x, δ)
            du = zeros(cachetype(x, δ), length(x))
            pbuf = Vector{cachetype(x, δ)}(p0)
            apply_perturb!(pbuf, δ, perturb_maps.p_map)
            nw(du, x, pbuf, t0; perturb=δ, perturb_maps)
            du
        end
        obsf = SII.observed(nw, outsym)
        g = function(x, δ)
            result = zeros(cachetype(x, δ), length(outsym))
            obsf(x, p0, t0, result; perturb=δ, perturb_maps)
            result
        end

        A = ForwardDiff.jacobian(x -> f(x, δ0), x0)
        B = ForwardDiff.jacobian(δ -> f(x0, δ), δ0)
        C = ForwardDiff.jacobian(x -> g(x, δ0), x0)
        D = ForwardDiff.jacobian(δ -> g(x0, δ), δ0)

        return NetworkDescriptorSystem(; M, A, B, C, D, sym,
            insym = single_in ? only(insym) : insym,
            outsym = single_out ? only(outsym) : outsym)
    end
end

"""
    reduce_dae(sys::NetworkDescriptorSystem) -> NetworkDescriptorSystem

Reduce a descriptor system `M ẋ = A x + B u`, `y = C x + D u` to an ODE system
by eliminating algebraic constraints. Returns a new `NetworkDescriptorSystem`
with `M = I` and reduced state dimension.

Partitions variables into differential (`M_ii = 1`) and algebraic (`M_ii = 0`).
The Jacobian `A` is partitioned as:

```
A = [f_x  f_y]    f_x = ∂f/∂x_d,  f_y = ∂f/∂x_a   (differential equations)
    [g_x  g_y]    g_x = ∂g/∂x_d,  g_y = ∂g/∂x_a   (algebraic constraints)
```

The algebraic constraint `0 = g_x*x_d + g_y*x_a + B_a*u` is solved for `x_a`
and substituted, giving the reduced system:

    A_s = f_x - f_y * g_y⁻¹ * g_x
    B_s = B_d - f_y * g_y⁻¹ * B_a
    C_s = C_d - C_a * g_y⁻¹ * g_x
    D_s = D   - C_a * g_y⁻¹ * B_a

For unconstrained systems (`M = I`), returns the system unchanged.

See: \"Power System Modelling and Scripting\", F. Milano, Chapter 7.2.
"""
function reduce_dae(sys::NetworkDescriptorSystem)
    M = sys.M
    A = sys.A
    if M isa LinearAlgebra.UniformScaling # treat as identity
        n = size(A, 1)
        d_idx = collect(1:n)
        return NetworkDescriptorSystem(;
            M = LinearAlgebra.I,
            A,
            B = sys.B,
            C = sys.C,
            D = sys.D,
            sym = sys.sym,
            insym = sys.insym,
            outsym = sys.outsym,
        )
    else
        diag_M = LinearAlgebra.diag(M)
        d_idx = findall(==(1), diag_M)
        c_idx = findall(==(0), diag_M)

        f_x = A[d_idx, d_idx]
        f_y = A[d_idx, c_idx]
        g_x = A[c_idx, d_idx]
        g_y = A[c_idx, c_idx]

        local gy_solve  # closure: x -> g_y \ x
        local A_s
        try
            gy_gx = g_y \ g_x
            gy_solve = x -> g_y \ x
            A_s = f_x - f_y * gy_gx
        catch e
            @warn "Matrix g_y is singular or ill-conditioned, falling back to pseudoinverse" exception=(e, catch_backtrace())
            gy_inv = LinearAlgebra.pinv(g_y)
            gy_solve = x -> gy_inv * x
            A_s = f_x - f_y * gy_inv * g_x
        end

        B_s = sys.B
        C_s = sys.C
        D_s = sys.D
        if !isnothing(sys.B)
            B_d = sys.B[d_idx, :]
            B_a = sys.B[c_idx, :]
            C_d = sys.C[:, d_idx]
            C_a = sys.C[:, c_idx]
            B_s = B_d - f_y * gy_solve(B_a)
            C_s = C_d - C_a * gy_solve(g_x)
            D_s = sys.D - C_a * gy_solve(B_a)
        end

        return NetworkDescriptorSystem(;
            M = LinearAlgebra.I,
            A = A_s,
            B = B_s,
            C = C_s,
            D = D_s,
            sym = isnothing(sys.sym) ? nothing : sys.sym[d_idx],
            insym = sys.insym,
            outsym = sys.outsym,
        )
    end
end

"""
    jacobian_eigenvals(s0::NWState; kwargs...)
    jacobian_eigenvals(sys::NetworkDescriptorSystem; eigvalf=LinearAlgebra.eigvals)

Compute the eigenvalues of the (reduced) Jacobian for linear stability analysis.

For unconstrained systems (`M = I`), returns eigenvalues of the full Jacobian.
For constrained systems (DAE), computes eigenvalues of the reduced Jacobian
after eliminating algebraic variables.

The `eigvalf` keyword allows passing a custom eigenvalue solver (e.g. an
autodifferentiable eigensolver).

See also: [`linearize_network`](@ref), [`participation_factors`](@ref)
"""
function jacobian_eigenvals(sys::NetworkDescriptorSystem; eigvalf=LinearAlgebra.eigvals)
    red = reduce_dae(sys)
    return eigvalf(red.A)
end

function jacobian_eigenvals(s0::NWState; kwargs...)
    return jacobian_eigenvals(linearize_network(s0); kwargs...)
end

"""
    is_linear_stable(s0; marginally_stable=false, atol=1e-14, kwargs...)

Check if the fixpoint `s0` is linearly stable by checking that all eigenvalues
of the (reduced) Jacobian have negative real parts.

Set `marginally_stable=true` to also accept eigenvalues with zero real part
(within tolerance `atol`).

See also: [`jacobian_eigenvals`](@ref), [`linearize_network`](@ref)
"""
function is_linear_stable(s0; marginally_stable=false, atol=1e-14, kwargs...)
    isfixpoint(s0) || error("The state s0 is not a fixpoint of the network.")
    λ = jacobian_eigenvals(s0; kwargs...)
    if marginally_stable
        return all(λi -> real(λi) < atol || isapprox(real(λi), 0.0; atol), λ)
    else
        return all(λi -> real(λi) < 0.0, λ)
    end
end

"""
    participation_factors(s0::NWState; kwargs...)
    participation_factors(sys::NetworkDescriptorSystem; normalize=true)


Compute participation factors showing which differential state variables
participate in each eigenmode of the (reduced) Jacobian.

Returns a named tuple:
- `eigenvalues`: eigenvalues of the reduced Jacobian
- `pfactors`: participation matrix where `pfactors[k,i]` is the participation
  of differential state `k` in mode `i`
- `sym`: symbolic names of the differential states

Participation factors are computed as `p_ki = |W[i,k] * V[k,i]|` where `V` are
right eigenvectors and `W = V⁻¹` are left eigenvectors. When `normalize=true`
(default), each mode's participation factors sum to 1.

See also: [`linearize_network`](@ref), [`jacobian_eigenvals`](@ref)
"""
function participation_factors(sys::NetworkDescriptorSystem; normalize=true)
    red = reduce_dae(sys)
    A_s = red.A

    eigenvalues, V = LinearAlgebra.eigen(A_s)
    W = LinearAlgebra.inv(V)

    n = size(A_s, 1)
    pfactors = zeros(n, n)
    for i in 1:n, k in 1:n
        pfactors[k, i] = abs(W[i, k] * V[k, i])
    end
    if normalize
        for i in 1:n
            s = sum(@view pfactors[:, i])
            s > 0 && (pfactors[:, i] ./= s)
        end
    end

    sym = isnothing(red.sym) ? nothing : red.sym
    return (; eigenvalues, pfactors, sym)
end

function participation_factors(s0::NWState; kwargs...)
    return participation_factors(linearize_network(s0); kwargs...)
end

"""
    show_participation_factors([io=stdout,] pf::NamedTuple; kwargs...)
    show_participation_factors(s0::NWState; normalize=true, kwargs...)

Display participation factors showing which states participate in each eigenmode.
Each mode is listed with its contributing states and their participation values.

When called with `s0`, internally calls [`participation_factors`](@ref) and then displays the
result. The `normalize` keyword is forwarded to `participation_factors`; remaining keywords control
display.

# Keyword Arguments
- `modes=nothing`: Filter which modes to display. Can be an integer (single mode),
  a range (`1:5`), a vector of indices, a predicate function on eigenvalues
  (e.g. `λ -> abs(imag(λ)) > 1`), or `nothing` (show all). Indices refer to
  the original eigenvalue ordering from `participation_factors`.
- `threshold=0.01`: Only show participation factors above this threshold (set to 0 to show all)
- `sigdigits=4`: Number of significant digits for displaying values
- `sortby=nothing`: Sort modes for display. `nothing` preserves `eigen`'s default
  lexicographic order (consistent with `participation_factors` and `eigenvalue_sensitivity`
  mode indices). Use `:eigenvalue` (real part, descending) or `:magnitude` for alternative orderings.

# Example
```julia
pf = participation_factors(s0)
show_participation_factors(pf)

# or directly:
show_participation_factors(s0; threshold=0.05)
```
"""
function show_participation_factors(io::IO, pf::NamedTuple;
                                    modes=nothing,
                                    threshold=0.10,
                                    sigdigits=4,
                                    sortby=nothing)
    (; eigenvalues, pfactors, sym) = pf

    # Select modes
    N = length(eigenvalues)
    selected = if isnothing(modes)
        collect(1:N)
    elseif modes isa Integer
        [modes]
    elseif modes isa Function
        findall(modes, eigenvalues)
    else
        collect(modes)
    end
    eigenvalues = eigenvalues[selected]
    pfactors = pfactors[:, selected]

    # Sort modes
    if sortby == :eigenvalue
        perm = sortperm(eigenvalues; by=real, rev=true)
    elseif sortby == :magnitude
        perm = sortperm(eigenvalues; by=abs, rev=true)
    else
        perm = 1:length(eigenvalues)
    end

    selected = selected[perm]
    eigenvalues = eigenvalues[perm]
    pfactors = pfactors[:, perm]

    n_modes = length(eigenvalues)
    n_states = size(pfactors, 1)

    # Build header
    header = "Mode & Real (Hz) & Imag (Hz) & Factor & State"

    # Build rows for each mode
    rows = String[]
    for i in 1:n_modes
        λ = eigenvalues[i]
        mode_label = string(selected[i])
        real_hz = real(λ) / (2π)
        imag_hz = imag(λ) / (2π)
        real_str = str_significant(real_hz; sigdigits, phantom_minus=true)
        imag_str = str_significant(imag_hz; sigdigits, phantom_minus=true)

        # Collect participating states with their values, sorted by participation
        participants = Tuple{Float64, String}[]
        for k in 1:n_states
            pf_val = pfactors[k, i]
            if abs(pf_val) >= threshold
                state_idx = sym === nothing ? "State($k)" : repr(sym[k])
                push!(participants, (pf_val, state_idx))
            end
        end
        sort!(participants; by=first, rev=true)

        # Create rows for this mode
        if isempty(participants)
            push!(rows, "$mode_label & $real_str & $imag_str & & -")
        else
            # First participation row includes eigenvalue
            pf_val, state_idx = participants[1]
            pf_str = str_significant(pf_val; sigdigits)
            push!(rows, "$mode_label & $real_str & $imag_str & $pf_str & $state_idx")
            # Subsequent participation rows have empty eigenvalue columns
            for (pf_val, state_idx) in participants[2:end]
                pf_str = str_significant(pf_val; sigdigits)
                push!(rows, " & & & $pf_str & $state_idx")
            end
        end
    end

    # Print table
    println(io, "Participation Factors:")
    if isempty(rows)
        println(io, "  (no modes)")
    else
        all_rows = [header; rows]
        aligned = align_strings(all_rows; padding=:right)
        println(io, "  ", aligned[1])
        println(io, "  ", "─"^maximum(textwidth.(aligned)))
        for line in aligned[2:end]
            println(io, "  ", line)
        end
    end
end
function show_participation_factors(pf::NamedTuple; kwargs...)
    show_participation_factors(stdout, pf; kwargs...)
end
function show_participation_factors(s0::NWState; kwargs...)
    show_participation_factors(stdout, s0; kwargs...)
end
function show_participation_factors(io, s0::NWState; normalize=true, kwargs...)
    pf = participation_factors(s0; normalize)
    show_participation_factors(io, pf; kwargs...)
end

"""
    eigenvalue_sensitivity(s0::NWState, mode_idx::Int; params=SII.parameter_symbols(extract_nw(s0)))

Compute the sensitivity of eigenvalue `mode_idx` to system parameters using
the classical eigenvalue sensitivity formula:

    ∂λ_i/∂p_k = wᵢᵀ · (∂Aₛ/∂p_k) · vᵢ

where `vᵢ` and `wᵢ` are the right and left eigenvectors of the reduced state
matrix `Aₛ` (with `W = V⁻¹`, so `wᵢᵀ·vᵢ = 1`). The Jacobian derivative
`∂Aₛ/∂p_k` is computed exactly via nested forward-mode automatic
differentiation (ForwardDiff).

Arguments:
- `mode_idx`: The indexes into eigenvalues in the default order of Julia's `eigen`
(lexicographic by `(real(λ), imag(λ))`). This is the same ordering as returned
by [`jacobian_eigenvals`](@ref) and [`participation_factors`](@ref).
- `params` (keyword): Symbolic parameter indices to compute sensitivities for. Defaults to
all parameters (`SII.parameter_symbols(s0)`). Pass a subset to reduce computation.

Returns a named tuple:
- `eigenvalue`: the selected eigenvalue (complex, rad/s)
- `mode_idx`: the mode index used
- `sensitivities`: `∂λ/∂p_k` for each selected parameter (complex vector)
- `scaled_sensitivities`: `p_k · ∂λ/∂p_k` for each selected parameter (eigenvalue shift per 100% parameter change)
- `param_syms`: symbolic parameter indices (matching the `params` argument)
- `param_vals`: parameter values at the operating point

See also: [`participation_factors`](@ref), [`jacobian_eigenvals`](@ref)
"""
function eigenvalue_sensitivity(s0::NWState, mode_idx::Int; params=SII.parameter_symbols(s0))
    nw = extract_nw(s0)
    # Baseline linearization and eigendecomposition
    A_s = reduce_dae(linearize_network(s0)).A

    eigenvalues, V = LinearAlgebra.eigen(A_s)
    W = LinearAlgebra.inv(V)

    n_modes = length(eigenvalues)
    1 <= mode_idx <= n_modes || throw(ArgumentError("mode_idx=$mode_idx out of range [1, $n_modes]"))

    λ = eigenvalues[mode_idx]
    v = V[:, mode_idx]       # right eigenvector
    w = W[mode_idx, :]       # left eigenvector (row of W = V⁻¹)
    # W = V⁻¹ implies wᵀ·v = 1 (no conjugation). Use transpose, not adjoint.
    # ∂λ/∂p = wᵀ · (∂A_s/∂p) · v  (denominator is 1 by construction)

    # Resolve parameter indices
    param_syms = if params isa SymbolicIndex && SII.symbolic_type(params) == SII.ScalarSymbolic()
        [params]
    else
        collect(params)
    end
    pidxs = SII.parameter_index(nw, param_syms)
    n_sel = length(param_syms)

    # Compute sensitivity for each selected parameter
    sensitivities = zeros(ComplexF64, n_sel)
    for (i, pidx) in enumerate(pidxs)
        dA_s = _reduced_jacobian_param_derivative(nw, s0, pidx)
        sensitivities[i] = transpose(w) * (dA_s * v)
    end

    param_vals = pflat(s0)[pidxs]
    scaled_sensitivities = map(1:n_sel) do i
        if isfinite(param_vals[i])
            param_vals[i] * sensitivities[i]
        else
            complex(NaN, NaN)
        end
    end

    return (; eigenvalue=λ, mode_idx, sensitivities, scaled_sensitivities,
              param_syms, param_vals)
end

function _reduced_jacobian_param_derivative(nw, s0, pidx)
    x0 = uflat(s0)
    p0 = pflat(s0)
    M = nw.mass_matrix
    ForwardDiff.derivative(0.0) do δp
        T = typeof(δp)
        pbuf = Vector{T}(p0)
        pbuf[pidx] += δp
        f(x) = begin
            du = zeros(promote_type(eltype(x), T), length(x))
            nw(du, x, pbuf, s0.t)
            du
        end
        J = ForwardDiff.jacobian(f, x0)
        reduce_dae(NetworkDescriptorSystem(; M, A=J)).A
    end
end

"""
    show_eigenvalue_sensitivity([io::IO, ]result::NamedTuple; threshold=0.01, sigdigits=3, sortby=:mag)
    show_eigenvalue_sensitivity([io::IO, ]s0::NWState, mode_idx; [params,] kwargs...)

Display eigenvalue sensitivity results in a formatted table, sorted by the specified metric
of the scaled sensitivity `p·∂λ/∂p`.

When called with `s0`, internally calls [`eigenvalue_sensitivity`](@ref) with `mode_idx` and
optional `params`, then displays the result. All keywords are forwarded to the display function.

# Keyword Arguments
- `threshold=0.01`: Only show entries where the chosen `sortby` metric exceeds this value (Hz or degrees for `:arg`)
- `sigdigits=4`: Number of significant digits
- `sortby=:mag`: Sort and filter by `:real`, `:imag`, `:mag`, or `:arg` of `p·∂λ/∂p`
"""
function show_eigenvalue_sensitivity(io::IO, result::NamedTuple;
                                     threshold=0.01, sigdigits=4, sortby=:mag)
    (; eigenvalue, mode_idx, sensitivities, scaled_sensitivities,
       param_syms, param_vals) = result

    λ = eigenvalue
    real_hz = real(λ) / (2π)
    imag_hz = imag(λ) / (2π)
    real_str = str_significant(real_hz; sigdigits, phantom_minus=true)
    imag_str = str_significant(imag_hz; sigdigits, phantom_minus=true)
    println(io, "Eigenvalue Sensitivity for mode $mode_idx: λ = $(real_str) ± $(abs(imag_hz) |> x -> str_significant(x; sigdigits))j Hz")

    # Collect and sort by chosen metric
    entries = Tuple{Float64, Int}[]
    for j in eachindex(scaled_sensitivities)
        s = scaled_sensitivities[j]
        sort_val = if sortby == :real
            real(s) / (2π)
        elseif sortby == :imag
            imag(s) / (2π)
        elseif sortby == :mag
            abs(s) / (2π)
        elseif sortby == :arg
            rad2deg(angle(s))
        elseif sortby == :realmag
            abs(real(s)) / (2π)
        elseif sortby == :imagmag
            abs(imag(s)) / (2π)
        else
            throw(ArgumentError("sortby must be :real, :imag, :mag, or :arg"))
        end

        if abs(sort_val) >= threshold
            push!(entries, (sort_val, j))
        end
    end
    sort!(entries; by=first, rev=true)

    if isempty(entries)
        unit = sortby == :arg ? "deg" : "Hz"
        println(io, "  (no sensitivities above threshold $threshold $unit)")
        return
    end

    header = "Parameter & Value & Real (Hz) & Imag (Hz) & |p·∂λ/∂p| (Hz) & ∠ (deg)"
    rows = String[]
    for (_, j) in entries
        sym_str = repr(param_syms[j])
        val_str = str_significant(param_vals[j]; sigdigits, phantom_minus=true)
        s = scaled_sensitivities[j]
        real_str = str_significant(real(s) / (2π); sigdigits, phantom_minus=true)
        imag_str = str_significant(imag(s) / (2π); sigdigits, phantom_minus=true)
        mag_str = str_significant(abs(s) / (2π); sigdigits)
        ang_str = str_significant(rad2deg(angle(s)); sigdigits, phantom_minus=true)
        push!(rows, "$sym_str & $val_str & $real_str & $imag_str & $mag_str & $ang_str")
    end

    all_rows = [header; rows]
    aligned = align_strings(all_rows; padding=:right)
    println(io, "  ", aligned[1])
    println(io, "  ", "─"^maximum(textwidth.(aligned)))
    for line in aligned[2:end]
        println(io, "  ", line)
    end
end
function show_eigenvalue_sensitivity(result::NamedTuple; kwargs...)
    show_eigenvalue_sensitivity(stdout, result; kwargs...)
end
function show_eigenvalue_sensitivity(io::IO, s0::NWState, mode_idx; params=SII.parameter_symbols(extract_nw(s0)), kwargs...)
    show_eigenvalue_sensitivity(io, eigenvalue_sensitivity(s0, mode_idx; params); kwargs...)
end
function show_eigenvalue_sensitivity(s0::NWState, mode_idx; kwargs...)
    show_eigenvalue_sensitivity(stdout, s0, mode_idx; kwargs...)
end

####
#### Internal: perturbation channel classification
####
function _classify_perturbation_channels(nw, perturbation_channels)
    δ0 = zeros(length(perturbation_channels))
    vo_map = Vector{NTuple{2, Int}}() # map outbuf idx -> δ idx
    vi_map = Vector{NTuple{2, Int}}() # map aggbuf idx -> δ idx
    eo_map = Vector{NTuple{2, Int}}() # map outbuf idx -> δ idx
    ei_map = Vector{NTuple{2, Int}}() # map gbuf idx -> δ idx
    p_map  = Vector{NTuple{2, Int}}() # map pflat idx -> δ idx

    for (i, sym) in enumerate(perturbation_channels)
        comp = getcomp(nw, sym)
        if SII.is_parameter(nw, sym)
            pidx = SII.parameter_index(nw, sym)
            push!(p_map, (pidx, i))
        elseif comp isa VertexModel
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
    perturb_maps = (; vo_map, vi_map, eo_map, ei_map, p_map)
    δ0, perturb_maps
end
