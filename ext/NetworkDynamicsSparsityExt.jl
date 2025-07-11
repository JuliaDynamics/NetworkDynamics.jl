module NetworkDynamicsSparsityExt

using SparseConnectivityTracer: TracerSparsityDetector, jacobian_sparsity
using NetworkDynamics: Network, NWState, uflat, pflat, resolvecompidx, ComponentModel,
                       EIndex, VIndex, EdgeModel, VertexModel,
                       StateMask, Symmetric, AntiSymmetric, Directed, Fiducial

"""
    get_jac_prototype(nw::Network; dense=false)

Compute the sparsity pattern of the Jacobian matrix for a NetworkDynamics network.

This function uses `SparseConnectivityTracer.jl` to detect the sparsity pattern of the Jacobian 
matrix of the network's dynamics function. The resulting sparsity pattern can be used to 
improve the performance of ODE solvers by providing structural information about the system.
The `dense` option is useful when certain components have complex sparsity patterns that
are difficult to detect automatically
The resulting prototype can significantly improve ODE solver performance for large networks

# Arguments
- `nw::Network`: The NetworkDynamics network for which to compute the Jacobian prototype
- `dense=false`: Controls which components should be treated as dense during sparsity detection:
  - `false`: Use actual component functions (default)
  - `true`: Replace all components with dense equivalents
  - `Vector{Union{VIndex, EIndex}}`: Replace only the specified vertex/edge components with dense equivalents

# Example Usage
```julia
nw = Network(...)
f_ode = ODEFunction(nw; jac_prototype=get_jac_prototype(nw))
prob = ODEProblem(f_ode, x0, (0.0, 1.0), p0)
sol = solve(prob, Rodas5P())
```
"""
function NetworkDynamics.get_jac_prototype(nw::Network; dense=false)
    if dense == true
        nw = replace_dense(nw, eachindex(nw.im.vertexm), eachindex(nw.im.edgem))
    elseif dense isa AbstractVector && !isempty(dense)
        vidxs = [resolvecompidx(nw, vidx) for vidx in dense if vidx isa VIndex]
        eidxs = [resolvecompidx(nw, eidx) for eidx in dense if eidx isa EIndex]
        nw = replace_dense(nw, vidxs, eidxs)
    elseif dense == false
        # do nothing
    else
        error("dense must be either true, false or a vector of VIndex and EIndex")
    end

    s0 = NWState(nw)
    x0 = uflat(s0)
    p0 = pflat(s0)
    fx = function(x)
        dx = similar(x)
        nw(dx, x, p0, 0.0)
        dx
    end

    detector = TracerSparsityDetector();
    jacobian_sparsity(fx, x0, detector)
end

function replace_dense(nw::Network, vidxs, eidxs)
    vertexm = copy(nw.im.vertexm)
    vertexm[vidxs] .= dense_equivalent.(nw.im.vertexm[vidxs])
    edgem = copy(nw.im.edgem)
    edgem[eidxs] .= dense_equivalent.(nw.im.edgem[eidxs])

    Network(nw; vertexm, edgem)
end

function dense_equivalent(cm::ComponentModel)
    ComponentModel = if cm isa VertexModel
        VertexModel
    elseif cm isa EdgeModel
        EdgeModel
    else
        error("ComponentModel must be either VertexModel or EdgeModel")
    end

    ComponentModel(cm; f=dense_f(cm), g=dense_g(cm), allow_output_sym_clash=true)
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

end # module
