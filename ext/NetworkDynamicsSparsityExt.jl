module NetworkDynamicsSparsityExt

using SparseConnectivityTracer: TracerSparsityDetector, jacobian_sparsity
using NetworkDynamics: Network, NWState, uflat, pflat, resolvecompidx, ComponentModel,
                       EIndex, VIndex, EdgeModel, VertexModel,
                       StateMask, Symmetric, AntiSymmetric, Directed, Fiducial

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
