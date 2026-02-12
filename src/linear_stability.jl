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

Without `in`/`out`, returns `(; M, A, state_syms)` containing the mass matrix,
state Jacobian, and symbolic state names.

With `in`/`out`, returns `(; M, A, B, C, D, G, state_syms, in_syms, out_syms)`
with full state-space matrices for the perturbation channels specified by `in`
and observations specified by `out`. `G(s)` is the transfer function
`G(s) = C * (M*s - A)⁻¹ * B + D`.

Perturbation channels (`in`) can be any component input or output, specified as
symbolic indices (e.g. `VIndex(1, :P)`, `EIndex(1, :i_r)`).

# Examples
```julia
# Simple linearization for eigenvalue analysis
sys = linearize_network(nw, s0)
eigvals = jacobian_eigenvals(sys)

# Full ABCD linearization with transfer function
sys = linearize_network(nw, s0; in=[VIndex(1, :u_r)], out=[VIndex(1, :i_r)])
sys.G(1im * 2π * 50)  # evaluate transfer function at 50 Hz
```

See also: [`jacobian_eigenvals`](@ref), [`participation_factors`](@ref)
"""
function linearize_network(nw, s0; in=nothing, out=nothing)
    x0 = uflat(s0)
    p0 = pflat(s0)
    t0 = s0.t
    M = nw.mass_matrix
    state_syms = SII.variable_symbols(nw)

    if isnothing(in) && isnothing(out)
        h!(dx, x) = nw(dx, x, p0, t0)
        A = ForwardDiff.jacobian(h!, similar(x0), x0)
        return (; M, A, state_syms)
    else
        isnothing(in) && throw(ArgumentError("Must specify `in` when `out` is given"))
        isnothing(out) && throw(ArgumentError("Must specify `out` when `in` is given"))

        in_syms = collect(in)
        out_syms = collect(out)
        δ0, perturb_maps = _classify_perturbation_channels(nw, in_syms)

        f = function(x, δ)
            du = zeros(promote_type(eltype(x), eltype(δ)), length(x))
            nw(du, x, p0, t0; perturb=δ, perturb_maps)
            du
        end
        obsf = SII.observed(nw, out_syms)
        g = function(x, δ)
            result = zeros(promote_type(eltype(x), eltype(δ)), length(out_syms))
            obsf(x, p0, t0, result; perturb=δ, perturb_maps)
            result
        end

        A = ForwardDiff.jacobian(x -> f(x, δ0), x0)
        B = ForwardDiff.jacobian(δ -> f(x0, δ), δ0)
        C = ForwardDiff.jacobian(x -> g(x, δ0), x0)
        D = ForwardDiff.jacobian(δ -> g(x0, δ), δ0)
        G(s) = C * ((M * s - A) \ B) + D

        return (; M, A, B, C, D, G, state_syms, in_syms, out_syms)
    end
end
linearize_network(s0::NWState, args...; kwargs...) = linearize_network(extract_nw(s0), s0, args...; kwargs...)

"""
    reduce_dae(M, A) -> (; A_s, d_idx, c_idx)

Reduce a differential-algebraic system `M * ẋ = A*x` to an ordinary differential
system by eliminating algebraic constraints.

Partitions variables into differential (`d_idx`, where `M_ii = 1`) and algebraic
(`c_idx`, where `M_ii = 0`). Let `x = [x_d; x_a]` where `x_d` are differential
and `x_a` algebraic variables. The Jacobian `A` is partitioned as:

```
A = [f_x  f_y]    f_x = ∂f/∂x_d,  f_y = ∂f/∂x_a   (differential equations)
    [g_x  g_y]    g_x = ∂g/∂x_d,  g_y = ∂g/∂x_a   (algebraic constraints)
```

The algebraic constraint `0 = g_x*x_d + g_y*x_a` is solved for `x_a` and
substituted into the differential equations, giving the reduced state matrix:

    A_s = f_x - f_y * g_y⁻¹ * g_x

For unconstrained systems (`M = I`), returns the original matrix unchanged.

See: \"Power System Modelling and Scripting\", F. Milano, Chapter 7.2.
"""
function reduce_dae(M, A)
    if M == LinearAlgebra.I
        n = size(A, 1)
        return (; A_s=A, d_idx=collect(1:n), c_idx=Int[])
    end

    diag_M = LinearAlgebra.diag(M)
    d_idx = findall(==(1), diag_M)
    c_idx = findall(==(0), diag_M)

    f_x = A[d_idx, d_idx]
    f_y = A[d_idx, c_idx]
    g_x = A[c_idx, d_idx]
    g_y = A[c_idx, c_idx]

    local A_s
    try
        A_s = f_x - f_y * (g_y \ g_x)
    catch e
        @warn "Matrix g_y is singular or ill-conditioned, falling back to pseudoinverse $g_y" exception=(e, catch_backtrace())
        gy_inv = LinearAlgebra.pinv(g_y)
        A_s = f_x - f_y * gy_inv * g_x
    end
    return (; A_s, d_idx, c_idx)
end

"""
    jacobian_eigenvals([nw::Network, ]s0::NWState; kwargs...)
    jacobian_eigenvals(sys::NamedTuple; eigvalf=LinearAlgebra.eigvals)

Compute the eigenvalues of the (reduced) Jacobian for linear stability analysis.

For unconstrained systems (`M = I`), returns eigenvalues of the full Jacobian.
For constrained systems (DAE), computes eigenvalues of the reduced Jacobian
after eliminating algebraic variables.

The `eigvalf` keyword allows passing a custom eigenvalue solver (e.g. an
autodifferentiable eigensolver).

See also: [`linearize_network`](@ref), [`participation_factors`](@ref)
"""
function jacobian_eigenvals(sys::NamedTuple; eigvalf=LinearAlgebra.eigvals)
    (; A_s) = reduce_dae(sys.M, sys.A)
    return eigvalf(A_s)
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
    participation_factors(sys::NamedTuple; normalize=true)

Compute participation factors showing which differential state variables
participate in each eigenmode of the (reduced) Jacobian.

Returns a named tuple:
- `eigenvalues`: eigenvalues of the reduced Jacobian
- `pfactors`: participation matrix where `pfactors[k,i]` is the participation
  of differential state `k` in mode `i`
- `state_syms`: symbolic names of the differential states

Participation factors are computed as `p_ki = |W[i,k] * V[k,i]|` where `V` are
right eigenvectors and `W = V⁻¹` are left eigenvectors. When `normalize=true`
(default), each mode's participation factors sum to 1.

See also: [`linearize_network`](@ref), [`jacobian_eigenvals`](@ref)
"""
function participation_factors(sys::NamedTuple; normalize=true)
    (; A_s, d_idx) = reduce_dae(sys.M, sys.A)

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

    state_syms = haskey(sys, :state_syms) ? sys.state_syms[d_idx] : nothing
    return (; eigenvalues, pfactors, state_syms)
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
- `threshold=0.01`: Only show participation factors above this threshold (set to 0 to show all)
- `sigdigits=3`: Number of significant digits for displaying values
- `sortby=:eigenvalue`: Sort modes by `:eigenvalue` (real part) or `:magnitude`

# Example
```julia
pf = participation_factors(nw, s0)
show_participation_factors(pf)
```
"""
function show_participation_factors(io::IO, pf::NamedTuple;
                                    threshold=0.10,
                                    sigdigits=3,
                                    sortby=:eigenvalue)
    (; eigenvalues, pfactors, state_syms) = pf

    # Sort modes
    if sortby == :eigenvalue
        perm = sortperm(eigenvalues; by=real, rev=true)
    elseif sortby == :magnitude
        perm = sortperm(eigenvalues; by=abs, rev=true)
    else
        perm = 1:length(eigenvalues)
    end

    eigenvalues = eigenvalues[perm]
    pfactors = pfactors[:, perm]

    n_modes = length(eigenvalues)
    n_states = size(pfactors, 1)

    # Build header
    header = "Real (Hz) & Imag (Hz) & Factor & State"

    # Build rows for each mode
    rows = String[]
    for i in 1:n_modes
        λ = eigenvalues[i]
        real_hz = real(λ) / (2π)
        imag_hz = imag(λ) / (2π)
        real_str = str_significant(real_hz; sigdigits, phantom_minus=true)
        imag_str = str_significant(imag_hz; sigdigits, phantom_minus=true)

        # Collect participating states with their values, sorted by participation
        participants = Tuple{Float64, String}[]
        for k in 1:n_states
            pf_val = pfactors[k, i]
            if abs(pf_val) >= threshold
                state_idx = state_syms === nothing ? "State($k)" : repr(state_syms[k])
                push!(participants, (pf_val, state_idx))
            end
        end
        sort!(participants; by=first, rev=true)

        # Create rows for this mode
        if isempty(participants)
            push!(rows, "$real_str & $imag_str & & -")
        else
            # First participation row includes eigenvalue
            pf_val, state_idx = participants[1]
            pf_str = str_significant(pf_val; sigdigits)
            push!(rows, "$real_str & $imag_str & $pf_str & $state_idx")
            # Subsequent participation rows have empty eigenvalue columns
            for (pf_val, state_idx) in participants[2:end]
                pf_str = str_significant(pf_val; sigdigits)
                push!(rows, " & & $pf_str & $state_idx")
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
