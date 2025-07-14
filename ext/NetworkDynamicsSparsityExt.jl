module NetworkDynamicsSparsityExt

using SparseConnectivityTracer: TracerSparsityDetector, jacobian_sparsity
using NetworkDynamics: NetworkDynamics, Network, NWState, uflat, pflat, resolvecompidx, ComponentModel,
                       EIndex, VIndex, EdgeModel, VertexModel, dim, pdim,
                       StateMask, Symmetric, AntiSymmetric, Directed, Fiducial
using MacroTools: @capture, postwalk
using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions, RuntimeGeneratedFunction, @RuntimeGeneratedFunction
using SparseArrays: sparse, findnz
using ForwardDiff: ForwardDiff
RuntimeGeneratedFunctions.init(@__MODULE__)

struct RemainingConditionalsException <: Exception end

"""
    get_jac_prototype(nw::Network; dense=false, remove_conditions=false)

Compute the sparsity pattern of the Jacobian matrix for a NetworkDynamics network.

This function uses `SparseConnectivityTracer.jl` to detect the sparsity pattern of the Jacobian 
matrix of the network's dynamics function. The resulting sparsity pattern can be used to 
improve the performance of ODE solvers by providing structural information about the system.
The `dense` option is useful when certain components have complex sparsity patterns that
are difficult to detect automatically. The `remove_conditions` option helps when conditional
statements in component functions interfere with sparsity detection.

# Arguments
- `nw::Network`: The NetworkDynamics network for which to compute the Jacobian prototype
- `dense=false`: Controls which components should be treated as dense during sparsity detection:
  - `false`: Use actual component functions (default)
  - `true`: Replace all components with dense equivalents
  - `Vector{Union{VIndex, EIndex}}`: Replace only the specified vertex/edge components with dense equivalents
- `remove_conditions=false`: Controls removal of conditional statements from component functions: this is only applicable
  to components defined via MTK. It essentially scans the function expression for if/else statements, deleting the condition
  and replacing the block by `truepath + falsepath`, which can help with sparsity detection.
  - `false`: Keep conditional statements as-is (default)
  - `true`: Remove conditionals from all components by converting `if-else` to additive form
  - `Vector{Union{VIndex, EIndex}}`: Remove conditionals only from specified vertex/edge components
- `check=true`: If `true`, the function checks the sparsity pattern against a forward-differentiated Jacobian as a sanity check.

# Returns
- A sparse matrix representing the sparsity pattern of the Jacobian matrix

# Example Usage
```julia
nw = Network(...)
jac_prototype = get_jac_prototype(nw) # get the sparsity pattern

# manually set define ODEFunction
f_ode = ODEFunction(nw; jac_prototype=jac_prototype)
prob = ODEProblem(f_ode, x0, (0.0, 1.0), p0)
sol = solve(prob, Rodas5P())

# ALTERNATIVE: use set_jac_prototype!
set_jac_prototype!(nw; jac_prototype) # attach pattern to network
prob = ODEProblem(nw, x0, (0.0, 1.0), p0) # uses jac prototype from network
sol = solve(prob, Rodas5P())
```

See also: [`set_jac_prototype!`](@ref)
"""
function NetworkDynamics.get_jac_prototype(nw::Network; dense=false, remove_conditions=false, check=true)
    nw_original = nw
    if dense == true
        nw = _replace_dense(nw, eachindex(nw.im.vertexm), eachindex(nw.im.edgem))
    elseif dense isa AbstractVector && !isempty(dense)
        vidxs = [resolvecompidx(nw, vidx) for vidx in dense if vidx isa VIndex]
        eidxs = [resolvecompidx(nw, eidx) for eidx in dense if eidx isa EIndex]
        nw = _replace_dense(nw, vidxs, eidxs)
    elseif dense == false
        # do nothing
    else
        error("`dense` must be either true, false or a vector of VIndex and EIndex")
    end

    if remove_conditions == true
        nw = filter_conditions(nw, eachindex(nw.im.vertexm), eachindex(nw.im.edgem))
    elseif remove_conditions isa AbstractVector && !isempty(remove_conditions)
        vidxs = [resolvecompidx(nw, vidx) for vidx in remove_conditions if vidx isa VIndex]
        eidxs = [resolvecompidx(nw, eidx) for eidx in remove_conditions if eidx isa EIndex]
        nw = filter_conditions(nw, vidxs, eidxs)
    elseif remove_conditions == false
        # do nothing
    else
        error("`remove_conditions` must be either true, false or a vector of VIndex and EIndex")
    end

    fx = function(xp)
        x = @views xp[1:dim(nw)]
        p = @views xp[dim(nw)+1:dim(nw)+pdim(nw)]
        dx = similar(x)
        nw(dx, x, p, 0.0)
        dx
    end

    s0 = NWState(nw)
    # to get a "global" pattern with respect to the parameters we need to determine
    # the jacobian for vcat(x0, p0)!
    detector = TracerSparsityDetector();
    jac = try
        jacobian_sparsity(fx, vcat(uflat(s0), pflat(s0)), detector)
    catch
        error("Automatic sparsity detection failed. Sometimes, this can be caused by a small \
               subset of components whose sparsity patterns are not detected correctly. You \
               can try the `remove_conditions` and `dense` options to ignore some of the \
               sub-sparsity patterns!")
    end
    ujac = jac[1:dim(nw), 1:dim(nw)] # slice off the upper part for the state variables

    if nw !== nw_original
        # map the symbols back to original ordering
        orig_symbols = NetworkDynamics.variable_symbols(nw_original)
        new_symbols = NetworkDynamics.variable_symbols(nw)
        index_map = map(ns -> findfirst(isequal(ns), orig_symbols), new_symbols)
    
        # Check that all symbols were found
        any(isnothing, index_map) && error("Could not map all symbols back to original ordering")

        # Extract I, J, V from the ujac sparse matrix
        I_new, J_new, V_new = findnz(ujac)

        # Map indices back to original symbol ordering
        I_orig = [index_map[i] for i in I_new]
        J_orig = [index_map[j] for j in J_new]

        # Construct the remapped sparse matrix
        n_orig = length(orig_symbols)
        ujac = sparse(I_orig, J_orig, V_new)
    end

    # compare with a forward diff jacobian to ensure sparsity pattern is correct
    if check
        s0_orig = NWState(nw_original)
        if all(!isnan, uflat(s0_orig)) && all(!isnan, pflat(s0_orig))
            fwjac = _jacobian_at_point(nw_original, s0_orig)
            _assert_conservative_pattern(fwjac, ujac)
        else
            s0_orig_ones = NWState(nw_original, ufill=1.0, pfill=1.0) # fill with ones
            fwjac = _jacobian_at_point(nw_original, s0_orig_ones)
            _assert_conservative_pattern(fwjac, ujac)
        end
    end
    return ujac
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


function _replace_dense(nw::Network, vidxs, eidxs)
    vertexm = copy(nw.im.vertexm)
    for idx in vidxs
        vertexm[idx] = dense_equivalent(nw.im.vertexm[idx])
    end
    edgem = copy(nw.im.edgem)
    for idx in eidxs
        edgem[idx] = dense_equivalent(nw.im.edgem[idx])
    end

    Network(nw; vertexm, edgem)
end

function dense_equivalent(cm::ComponentModel)
    comp_constructor(cm)(cm; f=dense_f(cm), g=dense_g(cm), allow_output_sym_clash=true)
end
function dense_f(cm)
    function(dx, args...)
         _dx = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
        for i in eachindex(dx)
            dx[i] = _dx
        end
        nothing
    end
end
function dense_g(cm::EdgeModel)
    if cm.g isa Union{
        Symmetric{<:StateMask},
        AntiSymmetric{<:StateMask},
        Directed{<:StateMask},
        Fiducial{<:StateMask, <:StateMask}
    }
        return cm.g
    end
    function(o1, o2, args...)
        _o = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
        for i in eachindex(o1)
            o1[i] = _o
        end
        for i in eachindex(o2)
            o2[i] = _o
        end
        nothing
    end
end
function dense_g(cm::VertexModel)
    if cm.g isa StateMask
        return cm.g
    end
    function(o1, args...)
        _o = mapreduce(arg -> isnothing(arg) ? 0.0 : sum(arg), +, args)
        for i in eachindex(o1)
            o1[i] = _o
        end
        nothing
    end
end

function filter_conditions(nw::Network, vidxs, eidxs)
    vertexm = copy(nw.im.vertexm)
    for idx in vidxs
        try
            vertexm[idx] = filter_conditions(nw.im.vertexm[idx])
        catch e
            if e isa RemainingConditionalsException
                @warn "NetworkDynamicsSparsityExt.filter_conditions: Remaining conditionals in vertex model $(idx)."
            else
                rethrow(e)
            end
        end
    end
    edgem = copy(nw.im.edgem)
    for idx in eidxs
        try
            edgem[idx] = filter_conditions(nw.im.edgem[idx])
        catch e
            if e isa RemainingConditionalsException
                @warn "NetworkDynamicsSparsityExt.filter_conditions: Remaining conditionals in edge model $(idx)."
            else
                rethrow(e)
            end
        end
    end
    Network(nw; vertexm, edgem)
end

function filter_conditions(cm::ComponentModel)
    f = filter_conditionals(cm.f)
    g = filter_conditionals(cm.g)
    if f!==cm.f || g!==cm.g
        comp_constructor(cm)(cm; f, g, allow_output_sym_clash=true)
    else
        cm
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
        if @capture(ex, if cond_; true_ex_; else; false_ex_; end)
            had_conditionals = true
            :(($true_ex) + ($false_ex))
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

function comp_constructor(cm)
    if cm isa VertexModel
        VertexModel
    elseif cm isa EdgeModel
        EdgeModel
    else
        error("ComponentModel must be either VertexModel or EdgeModel")
    end
end

end # module
