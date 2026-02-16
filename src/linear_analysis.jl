"""
    NetworkDescriptorSystem

Descriptor system representation

    M ẋ = A x + B u
      y = C x + D u

of a linearized network. Returned by [`linearize_network`](@ref).

The struct is callable as a transfer function: `sys(s)` computes
`C (M s - A)⁻¹ B + D`. For SISO systems (scalar `in_syms`/`out_syms`) the
result is a scalar.

Acces the properties using `sys.A`, `.B`, `.C`, `.D`, `.M` to obtain matrices.
Use `sym(⋅)`, `insym(⋅)` and `outsym(⋅)` to obtain state names.

See also: [`linearize_network`](@ref), [`reduce_dae`](@ref)
"""
@kwdef struct NetworkDescriptorSystem{MT,AT,BT,CT,DT,IST,OST}
    M::MT = LinearAlgebra.I
    A::AT
    B::BT = nothing
    C::CT = nothing
    D::DT = nothing
    sym::Vector
    insym::IST = nothing
    outsym::OST = nothing
    function NetworkDescriptorSystem(M, A, B, C, D, sym, insym, outsym)
        dim = size(A, 1)
        if length(sym) != dim
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

function Base.show(io::IO, ::MIME"text/plain", sys::NetworkDescriptorSystem)
    compact = get(io, :compact, false)::Bool
    printstyled(io, "NetworkDescriptorSystem"; bold=true)
    if sys.B === nothing
        if compact
            print(io, "(dim=$(dim(sys)))")
        else
            print(io, " with\n")
            print(io, "  ├─ $(dim(sys)) states\n")
            print(io, "  └─ without inputs/outputs\n")
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
    isfixpoint([nw::Network, ]s0::NWState; tol=1e-10)

Check if the state `s0` is a fixpoint of the network `nw` by
calculating the the RHS and check that every entry is within the
given tolerance `tol`.
"""
function isfixpoint(nw::Network, s0::NWState; tol=1e-10)
    u0 = uflat(s0)
    p0 = pflat(s0)
    du = zeros(eltype(u0), length(u0))
    nw(du, u0, p0, s0.t)
    return all(abs.(du) .< tol)
end
isfixpoint(s0::NWState, args...; kwargs...) = isfixpoint(extract_nw(s0), s0, args...; kwargs...)

"""
    linearize_network([nw::Network, ]s0::NWState)
    linearize_network([nw::Network, ]s0::NWState; in, out)

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
sys = linearize_network(nw, s0)
eigvals = jacobian_eigenvals(sys)

# Full ABCD linearization with transfer function
sys = linearize_network(nw, s0; in=[VIndex(1, :u_r)], out=[VIndex(1, :i_r)])
sys(im * 2π * 50)  # evaluate transfer function at 50 Hz
```

See also: [`jacobian_eigenvals`](@ref), [`participation_factors`](@ref)
"""
function linearize_network(nw, s0; in=nothing, out=nothing)
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
            du = zeros(promote_type(eltype(x), eltype(δ)), length(x))
            pbuf = Vector{promote_type(eltype(x), eltype(δ))}(p0)
            apply_perturb!(pbuf, δ, perturb_maps.p_map)
            nw(du, x, pbuf, t0; perturb=δ, perturb_maps)
            du
        end
        obsf = SII.observed(nw, outsym)
        g = function(x, δ)
            result = zeros(promote_type(eltype(x), eltype(δ)), length(outsym))
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
linearize_network(s0::NWState, args...; kwargs...) = linearize_network(extract_nw(s0), s0, args...; kwargs...)

"""
    reduce_dae(sys::NetworkDescriptorSystem) -> NetworkDescriptorSystem

Reduce a descriptor system `M ẋ = A x` to an ODE system by eliminating algebraic
constraints. Returns a new `NetworkDescriptorSystem` with `M = I` and reduced
state dimension.

Partitions variables into differential (`M_ii = 1`) and algebraic (`M_ii = 0`).
The Jacobian `A` is partitioned as:

```
A = [f_x  f_y]    f_x = ∂f/∂x_d,  f_y = ∂f/∂x_a   (differential equations)
    [g_x  g_y]    g_x = ∂g/∂x_d,  g_y = ∂g/∂x_a   (algebraic constraints)
```

The algebraic constraint `0 = g_x*x_d + g_y*x_a` is solved for `x_a` and
substituted into the differential equations, giving the reduced state matrix:

    A_s = f_x - f_y * g_y⁻¹ * g_x

For unconstrained systems (`M = I`), returns the system with unchanged `A`.

See: \"Power System Modelling and Scripting\", F. Milano, Chapter 7.2.
"""
function reduce_dae(sys::NetworkDescriptorSystem)
    M = sys.M
    A = sys.A
    if M == LinearAlgebra.I
        n = size(A, 1)
        d_idx = collect(1:n)
    else
        diag_M = LinearAlgebra.diag(M)
        d_idx = findall(==(1), diag_M)
        c_idx = findall(==(0), diag_M)

        f_x = A[d_idx, d_idx]
        f_y = A[d_idx, c_idx]
        g_x = A[c_idx, d_idx]
        g_y = A[c_idx, c_idx]

        try
            A = f_x - f_y * (g_y \ g_x)
        catch e
            @warn "Matrix g_y is singular or ill-conditioned, falling back to pseudoinverse $g_y" exception=(e, catch_backtrace())
            gy_inv = LinearAlgebra.pinv(g_y)
            A = f_x - f_y * gy_inv * g_x
        end
    end
    return NetworkDescriptorSystem(;
        M = LinearAlgebra.I,
        A,
        B = sys.B,
        C = sys.C,
        D = sys.D,
        sym = sys.sym[d_idx],
        insym = sys.insym,
        outsym = sys.outsym,
    )
end

"""
    jacobian_eigenvals([nw::Network, ]s0::NWState; kwargs...)
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

function jacobian_eigenvals(nw::Network, s0; kwargs...)
    return jacobian_eigenvals(linearize_network(nw, s0); kwargs...)
end
jacobian_eigenvals(s0::NWState, args...; kwargs...) = jacobian_eigenvals(extract_nw(s0), s0, args...; kwargs...)

"""
    is_linear_stable([nw::Network, ]s0; marginally_stable=false, atol=1e-14, kwargs...)

Check if the fixpoint `s0` is linearly stable by checking that all eigenvalues
of the (reduced) Jacobian have negative real parts.

Set `marginally_stable=true` to also accept eigenvalues with zero real part
(within tolerance `atol`).

See also: [`jacobian_eigenvals`](@ref), [`linearize_network`](@ref)
"""
function is_linear_stable(nw::Network, s0; marginally_stable=false, atol=1e-14, kwargs...)
    isfixpoint(nw, s0) || error("The state s0 is not a fixpoint of the network nw.")
    λ = jacobian_eigenvals(nw, s0; kwargs...)
    if marginally_stable
        return all(λi -> real(λi) < atol || isapprox(real(λi), 0.0; atol), λ)
    else
        return all(λi -> real(λi) < 0.0, λ)
    end
end
is_linear_stable(s0::NWState, args...; kwargs...) = is_linear_stable(extract_nw(s0), s0, args...; kwargs...)

"""
    participation_factors([nw::Network, ]s0::NWState; kwargs...)
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

    sym = isempty(red.sym) ? nothing : red.sym
    return (; eigenvalues, pfactors, sym)
end

function participation_factors(nw::Network, s0; kwargs...)
    return participation_factors(linearize_network(nw, s0); kwargs...)
end
participation_factors(s0::NWState, args...; kwargs...) = participation_factors(extract_nw(s0), s0, args...; kwargs...)

"""
    show_participation_factors(io::IO, pf::NamedTuple; kwargs...)
    show_participation_factors(pf::NamedTuple; kwargs...)

Display participation factors showing which states participate in each eigenmode.
Each mode is listed with its contributing states and their participation values.

# Keyword Arguments
- `modes=nothing`: Filter which modes to display. Can be an integer (single mode),
  a range (`1:5`), a vector of indices, a predicate function on eigenvalues
  (e.g. `λ -> abs(imag(λ)) > 1`), or `nothing` (show all). Indices refer to
  the original eigenvalue ordering from `participation_factors`.
- `threshold=0.01`: Only show participation factors above this threshold (set to 0 to show all)
- `sigdigits=3`: Number of significant digits for displaying values
- `sortby=nothing`: Sort modes for display. `nothing` preserves `eigen`'s default
  lexicographic order (consistent with `participation_factors` and `eigenvalue_sensitivity`
  mode indices). Use `:eigenvalue` (real part, descending) or `:magnitude` for alternative orderings.

# Example
```julia
pf = participation_factors(nw, s0)
show_participation_factors(pf)
```
"""
function show_participation_factors(io::IO, pf::NamedTuple;
                                    modes=nothing,
                                    threshold=0.10,
                                    sigdigits=3,
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
function show_participation_factors(s0::NWState, args...; kwargs...)
    show_participation_factors(extract_nw(s0), s0, args...; kwargs...)
end
function show_participation_factors(nw::Network, s0; kwargs...)
    pf = participation_factors(nw, s0)
    show_participation_factors(pf; kwargs...)
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
