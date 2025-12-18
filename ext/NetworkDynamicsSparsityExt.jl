module NetworkDynamicsSparsityExt

using SparseConnectivityTracer: TracerSparsityDetector, jacobian_sparsity
using NetworkDynamics: NetworkDynamics, Network, NWState, uflat, pflat,
                       EdgeModel, VertexModel, dim, pdim,
                       Symmetric, AntiSymmetric, Directed, Fiducial,
                       executionstyle, SequentialAggregator, Aggregator, KAAggregator,
                       ComponentBatch, indim, outdim, extdim, fftype
using MacroTools: @capture, postwalk
using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions, RuntimeGeneratedFunction, @RuntimeGeneratedFunction
using ForwardDiff: ForwardDiff
using Accessors: @set
RuntimeGeneratedFunctions.init(@__MODULE__)

struct RemainingConditionalsException <: Exception end

"""
    get_jac_prototype(nw::Network; check=:auto, verbose=true)

Compute the sparsity pattern of the Jacobian matrix for a NetworkDynamics network.

This function uses `SparseConnectivityTracer.jl` (SCT) to detect the sparsity pattern of the Jacobian
matrix of the network's dynamics function. The resulting sparsity pattern can be used to 
improve the performance of ODE solvers by providing structural information about the system.

On a per-batch basis (i.e. once per unique component), the function will attempt to get the
sparsity pattern. Sometimes, certain component functions may not be compatible with SCT,
in that case the function will attempt to
- first try to replace if/else in RuntimeGeneratedFunctions with ifelse-function
- if it still does not work, replace the component function with a dense equivalent.

The `check` argument controls whether the detected sparsity pattern is verified against
the forward-diff jacobian. This check will be disabled with a warning for large networks.

# Returns
- A sparse matrix representing the sparsity pattern of the Jacobian matrix

See also [`NetworkDynamics.set_jac_prototype!`](@ref).
"""
function NetworkDynamics.get_jac_prototype(
    nw_original::Network;
    check=:auto, make_compatible=true, verbose=true,
    dense=nothing, remove_conditions=nothing, # deprecated
)
    if check == :auto
        if dim(nw_original) > 10_000
            check = false
            verbose && @warn "Skipping sparsity pattern check for large network (dim=$(dim(nw_original))). Force using `check=true`"
        else
            check = true
        end
    end

    if !isnothing(dense)
        @warn "The `dense` keyword argument for `get_jac_prototype` is deprecated. The algorithm will automatically fall back to this option on a per-batch basis."
    end
    if !isnothing(remove_conditions)
        @warn "The `remove_conditions` keyword argument for `get_jac_prototype` is deprecated. The algorithm will automatically fall back to this option on a per-batch basis."
    end

    nw_compat = make_compatible ? _to_compatible_execution_scheme(nw_original) : nw_original
    vbatches = get_compatible_batch.(nw_compat.vertexbatches; verbose);
    ebatches = get_compatible_batch.(nw_compat.layer.edgebatches; verbose);

    # replace both vertex and edge batches with compatible versions
    _nw = @set nw_compat.vertexbatches = vbatches
    nw = @set _nw.layer.edgebatches = ebatches

    s0 = NWState(nw)
    wrap = let nw=nw, p=pflat(s0)
        function(dx, x)
            nw(dx, x, p, NaN)
            nothing
        end
    end

    detector = TracerSparsityDetector()
    jac = jacobian_sparsity(wrap, similar(uflat(s0)), uflat(s0), detector)

    # compare with a forward diff jacobian to ensure sparsity pattern is correct
    if check
        s0_orig = NWState(nw_compat)
        if all(!isnan, uflat(s0_orig)) && all(!isnan, pflat(s0_orig))
            fwjac = _jacobian_at_point(nw_compat, s0_orig)
            _assert_conservative_pattern(fwjac, jac)
        else
            s0_orig_ones = NWState(nw_compat, ufill=1.0, pfill=1.0) # fill with ones
            fwjac = _jacobian_at_point(nw_compat, s0_orig_ones)
            _assert_conservative_pattern(fwjac, jac)
        end
    end
    return jac
end

function get_compatible_batch(batch::ComponentBatch{T}; verbose) where T
    (;fworks, gworks) = _test_SCT_compat(batch)
    fworks && gworks && return batch

    if T == VertexModel
        type = "vertex batch"
    elseif T == EdgeModel
        type = "edge batch"
    else
        error("ComponentBatch type must be either VertexModel or EdgeModel")
    end

    verbose && println("There was a problem in $type containing $(length((batch.indices))) components:")

    # stage 1, try to replace conditionals
    retest = false
    if !fworks && batch.compf isa RuntimeGeneratedFunction
        verbose && println(" - Try to filter conditionals in RuntimeGeneratedFunction compf")
        local newf
        try
            newf = filter_conditionals(batch.compf)
        catch e
            verbose && println(" - Failed to remove conditionals from compf: $e")
            newf = nothing
        end
        if !isnothing(newf)
            batch = @set batch.compf = newf
            retest = true
        end
    end
    if !gworks
        local newg
        try
            newg = filter_conditionals(batch.compg) # only does something for RGF wrappers
        catch e
            newg = nothing
        end
        if newg !== batch.compg #
            verbose && println(" - Try to filter conditionals in RuntimeGeneratedFunctions compg")
            if isnothing(newg)
                verbose && println(" - Failed to remove conditionals from compg: $e")
            else
                batch = @set batch.compg = newg
                retest = true
            end
        end
    end

    if retest # only retest if something changed
        (;fworks, gworks) = _test_SCT_compat(batch; error=false)
    end
    fworks && gworks && return batch

    # stage 2, replace with dense equivalents
    if !fworks
        verbose && println(" - Replacing compf with dense equivalent")
        batch = @set batch.compf = _generic_dense_f
    end
    if !gworks
        verbose && println(" - Replacing compg with dense equivalent")
        densg = T == VertexModel ? _generic_dense_vertexg : _generic_dense_edgeg
        batch = @set batch.compg = densg
    end

    (;fworks, gworks) = _test_SCT_compat(batch; error=true)
    fworks && gworks && return batch

    error("Could not make $type compatible with SparseConnectivityTracer.jl")
end
function _test_SCT_compat(batch; error=false)
    fworks = false
    gworks = false

    xdim = dim(batch)

    urange = 1:Int(dim(batch))
    inranges = UnitRange{Int}[]
    last_inrange_idx = dim(batch)
    allindims = extdim(batch) > 0 ? (indim(batch)..., extdim(batch)) : (indim(batch)...,)
    for idim in allindims
        push!(inranges, last_inrange_idx+1:last_inrange_idx+idim)
        last_inrange_idx += idim
    end
    prange = last_inrange_idx .+ (1 : pdim(batch))

    detector = TracerSparsityDetector()
    flatin = zeros(last_inrange_idx + pdim(batch) + 1)

    if isnothing(batch.compf)
        fworks=true
    else
        batchf = (du, allinputs) -> begin
            u = view(allinputs, urange)
            ins = map(i -> view(allinputs, inranges[i]), 1:length(inranges))
            p = view(allinputs, prange)
            NetworkDynamics.apply_compf(batch.compf, du, u, ins, p, allinputs[end])
            nothing
        end

        dx = zeros(Int(xdim))
        try
            jac = jacobian_sparsity(batchf, dx, flatin, detector)
            fworks = true
        catch e
            # @error "while f" e
            error && rethrow(e)
        end
    end

    outranges = UnitRange{Int}[]
    last_outrange_idx = 0
    for odim in outdim(batch)
        push!(outranges, last_outrange_idx+1:last_outrange_idx+odim)
        last_outrange_idx += odim
    end

    batchg = (allgoutputs, allinputs) -> begin
        outs = map(i -> view(allgoutputs, outranges[i]), 1:length(outranges))
        u = view(allinputs, urange)
        ins = map(i -> view(allinputs, inranges[i]), 1:length(inranges))
        p = view(allinputs, prange)
        NetworkDynamics.apply_compg(fftype(batch), batch.compg, outs, u, ins, p, allinputs[end])
        nothing
    end

    flatout = zeros(last_outrange_idx)
    try
        jac = jacobian_sparsity(batchg, flatout, flatin, detector)
        gworks = true
    catch e
        # @error "while g" e
        error && rethrow(e)
    end

    return (; fworks, gworks)
end

function _jacobian_at_point(nw, s0)
    ForwardDiff.jacobian(
        (du, u) -> nw(du, u, pflat(s0), NaN),
        zeros(dim(nw)),
        uflat(s0)
    )
end

function _assert_conservative_pattern(refjac, pattern)
    @assert size(refjac) == size(pattern) "Reference jacobian and template jacobian must have the same size!"
    for i in eachindex(refjac)
        if !iszero(refjac[i]) && iszero(pattern[i])
            error("Sparsity pattern mismatch! ForwardDiff returned nonzero entry where pattern is zero!")
        end
    end
end

function filter_conditionals(rgf::RuntimeGeneratedFunction)
    ex = RuntimeGeneratedFunctions.get_expression(rgf)
    newex = filter_conditionals_expr(ex)
    if newex !== ex
        return @RuntimeGeneratedFunction(newex)
    else
        return rgf
    end
end
filter_conditionals(x) = x

function filter_conditionals(w::AntiSymmetric)
    g = filter_conditionals(w.g)
    g !== w.g ? AntiSymmetric(g) : w
end
function filter_conditionals(w::Symmetric)
    g = filter_conditionals(w.g)
    g !== w.g ? Symmetric(g) : w
end
function filter_conditionals(w::Directed)
    g = filter_conditionals(w.g)
    g !== w.g ? Directed(g) : w
end
function filter_conditionals(w::Fiducial)
    dst = filter_conditionals(w.dst)
    src = filter_conditionals(w.src)
    if dst !== w.dst || src !== w.src
        Fiducial(dst, src)
    else
        w
    end
end

function filter_conditionals_expr(expr)
    had_conditionals = false
    newex = postwalk(expr) do ex
        if @capture(ex, lhs_ = if cond_; true_ex_; else; false_ex_; end)
            had_conditionals = true
            :($lhs = ifelse($cond, $true_ex, $false_ex))
        else
            ex
        end
    end

    remaining_ifs = false
    postwalk(newex) do ex
        if ex isa Expr && (ex.head == :if || ex.head == :elseif)
            remaining_ifs = true
        end
        nothing
    end
    # return newex

    if !had_conditionals && !remaining_ifs
        return expr
    elseif had_conditionals && !remaining_ifs
        return newex
    else
        throw(RemainingConditionalsException())
    end
end

function _to_compatible_execution_scheme(nw_original)
    newex = _compatible_exstyle(nw_original)
    if newex !== executionstyle(nw_original) || _needs_new_aggregator(nw_original)
        newagg = _new_aggregator(nw_original.layer.aggregator)
        return Network(nw_original; execution=newex, aggregator=newagg)
    else
        return nw_original
    end
end
# all ex styles are compatible
_compatible_exstyle(nw) = executionstyle(nw)

# KAAggregator is known to be incompatible
_needs_new_aggregator(nw) = nw.layer.aggregator isa KAAggregator
_new_aggregator(x::KAAggregator) = SequentialAggregator(x.f)
_new_aggregator(x::Aggregator) = NetworkDynamics.get_aggr_constructor(x)


function _generic_dense_f(dx, args...)
        _dx = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
    for i in eachindex(dx)
        dx[i] = _dx
    end
    nothing
end
function _generic_dense_edgeg(o1, o2, args...)
    _o = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
    for i in eachindex(o1)
        o1[i] = _o
    end
    for i in eachindex(o2)
        o2[i] = _o
    end
    nothing
end
function _generic_dense_vertexg(o1, args...)
    _o = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
    for i in eachindex(o1)
        o1[i] = _o
    end
    nothing
end

end # module
