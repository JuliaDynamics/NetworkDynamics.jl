abstract type ComponentCallback{C,A} end

struct ComponentCondition{C,DIM,PDIM}
    f::C
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
end

struct ComponentAffect{A,DIM,PDIM}
    f::A
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
end

# ComponentCondition([:u1, :u2], [:p1, :p2]) do out, u, p, t
#     out[1] = u[1] + p[1]
#     out[2] = u[2] + p[2]
# end

# ComponentAffect([], [:p1]) do u, p, idx, ctx
#     p[:p1] = 0
# end


struct ContinousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback{C,A}
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    kwargs::KW
end

struct VectorContinousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback{C,A}
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    len::Int
    kwargs::KW
end

struct CallbackBatch{T,C,A}
    nw::Network
    components::Vector{SymbolicIndex}
    callbacks::Vector{T}
    sublen::Int # length of each callback
    condition::C
    affect::A
end
function CallbackBatch(nw, components, callbacks)
    if !isconcretetype(eltype(componens))
        components = Vector{typeof(first(components))}(components)
    end
    if !isconcretetype(eltype(componens))
        callbacks = Vector{typeof(first(callbacks))}(callbacks)
    end
    sublen = eltype(components) isa ContinousComponentCallback ? 1 : first(components).len
    condition = first(callbacks).condition
    affect = first(callbacks).affect
    CallbackBatch(nw, components, callbacks, sublen, condition, affect)
end

condition_dim(cbb)  = first(cbb.callbacks).condition.sym  |> length
condition_pdim(cbb) = first(cbb.callbacks).condition.psym |> length
affect_dim(cbb)  = first(cbb.callbacks).affect.sym  |> length
affect_pdim(cbb) = first(cbb.callbacks).affect.psym |> length

condition_urange(cbb, i) = (1 + (i-1)*condition_dim(cbb)  + 1) : i*condition_dim(cbb)
condition_prange(cbb, i) = (1 + (i-1)*condition_pdim(cbb) + 1) : i*condition_pdim(cbb)
affect_urange(cbb, i) = (1 + (i-1)*affect_dim(cbb)  + 1) : i*affect_dim(cbb)
affect_prange(cbb, i) = (1 + (i-1)*affect_pdim(cbb) + 1) : i*affect_pdim(cbb)

condition_outrange(cbb, i) = (1 + (i-1)*cbb.sublen + 1) : i*cbb.sublen

cbidx_from_outidx(cbb, outidx) = div(cbb.sublen, outidx) + 1

function collect_c_or_a_indices(cbb, c_or_a, u_or_p)
    for component, cb in zip(cbb.components, cbb.callbacks)
        syms = getproperty(getproperty(cb, c_or_a), u_or_p)
        sidx = collect(typeof(component)(component.compidx, syms))
        append!(sidxs, sidx)
    end
end

function collect_callbackbatches(nw)
    component = SymbolicIndex[]
    callbacks = ComponentCallback[]
    for (i, v) in pairs(nw.im.vertexm)
        has_callback(v) || continue
        for cb in get_callback(v)
            push!(components, VIndex(i, nothing))
            push!(callbacks, cb)
        end
    end
    for (i, v) in pairs(nw.im.edgem)
        has_callback(v) || continue
        for cb in get_callback(v)
            push!(components, EIndex(i, nothing))
            push!(callbacks, cb)
        end
    end

    idx_per_type = _find_identical(components, 1:length(components))
    batches = CallbackBatch[]
    for typeidx in idx_per_type
        batchcomps = components[typeidx]
        batchcbs = callbacks[typeidx]
        cbs = CallbackBatch(nw, batchcomps, batchcbs)
        push!(batches, cbs)
    end
    return cbs
end

function batchequal(a::ContinousComponentCallback, b::ContinousComponentCallback)
    batchequal(a.condition, b.condition) || return false
    batchequal(a.affect, b.affect)       || return false
    batchequal(a.kwargs, b.kwargs)       || return false
    return true
end
function batchequal(a::VectorContinousComponentCallback, b::VectorContinousComponentCallback)
    batchequal(a.condition, b.condition) || return false
    batchequal(a.affect, b.affect)       || return false
    batchequal(a.kwargs, b.kwargs)       || return false
    batchequal(a.len, b.len)             || return false
    return true
end
function batchequal(a::T, b::T) where {T <: Union{ComponentCondition, ComponentAffect}}
    typeof(a) == typeof(b)
end
function batchequal(a::NamedTuple, b::NamedTuple)
    length(a) == length(b) || return false
    for (k, v) in a
        haskey(b, k) || return false
        v == b[k] || return false
    end
    return true
end


function c_symbolic_view_p(cbb::CallbackBatch, p, i)
    _x = view(p, condition_prange(cbb, i))
    SymbolicView(_x, cbb.callbacks[i].condition.psym)
end
function c_symbolic_view_u(cbb::CallbackBatch, u, i)
    _x = view(u, condition_urange(cbb, i))
    SymbolicView(_x, cbb.callbacks[i].condition.sym)
end
function a_symbolic_view_p(cbb::CallbackBatch, p, i)
    _x = view(p, affect_prange(cbb, i))
    SymbolicView(_x, cbb.callbacks[i].affect.psym)
end
function a_symbolic_view_u(cbb::CallbackBatch, u, i)
    _x = view(u, affect_urange(cbb, i))
    SymbolicView(_x, cbb.callbacks[i].affect.sym)
end

function batch_condition(cbbatch)
    ucache = DiffCache(zeros(dim(cbatch)), 12)
    usymidxs = collect_c_or_a_indices(cbb, :condition, :sym)
    psymidxs = collect_c_or_a_indices(cbb, :condition, :psym)
    obsf = SII.observed(cbbatch.nw, usymidxs)
    pidxs = SII.parameter_index(cbbatch.nw, psymidxs)

    (out, u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us)
        ps = @view integrator.p[pidxs]

        for i in 1:length(cbbatch)
            _u = c_symbolic_view_u(cbbatch, us, i)
            _p = c_symbolic_view_p(cbbatch, ps, i)

            if cbtype(cbbatch) <: VectorContinousComponentCallback
                _out = batch.condition(_u, _p, t)
            elseif cbtype(cbbatch) <: ContinousComponentCallback
                @views _out = out[condition_outrange(cbbatch, i)]
                cbbatch.condition(_out, _u, _p, t)
            end
        end
        nothing
    end
end
function batch_affect(cbb)
    usymidxs = collect_c_or_a_indices(cbb, :affect, :sym)
    psymidxs = collect_c_or_a_indices(cbb, :affect, :psym)
    uidxs = SII.variable_index(cbb.nw, usymidx)
    pidxs = SII.parameter_index(cbb.nw, psymidx)

    (integrator, outidx) -> begin
        us = @view integrator.u[uidxs]
        ps = @view integrator.p[pidxs]

        i = cbidx_from_outidx(cbb, outidx)
        _u = a_symbolic_view_u(cbb, us, i)
        _p = a_symbolic_view_p(cbb, ps, i)

        uhash = hash(_u.v)
        phash = hash(_p.v)
        cbb.affect(_u, _p, integrator.t, get_ctx(cbb, i))
        pchanged = hash(_p.v) != phash
        uchanged = hash(_u.v) != uhash

        (pchanged || uchanged) && auto_dt_reset!(integrator)
        pchanged && save_parameters!(integrator)
    end
end

get_ctx(cbb, i::Int) = get_ctx(cbb, cbb.components[i])
function get_ctx(cbb, sym::VIndex)
    idx = sym.compidx
    (; model=cbb.nw.im.vertexm[idx], vidx=idx)
end
function get_ctx(cbb, sym::EIndex)
    idx = sym.compidx
    edge = cbb.nw.im.edgevec[idx]
    (; model=cbb.nw.im.edgem[idx], eidx=idx, src=edge.src, dst=edge.dst)
end

struct SymbolicView{N,VT}
    v::VT
    syms::NTuple{N,Symbol}
end
Base.IteratorSize(::Type{SymbolicView}) = IteratorSize(x.v)
Base.IteratorEltype(::Type{SymbolicView}) = IteratorEltype(x.v)
Base.eltype(::Type{SymbolicView}) = eltype(x.v)
Base.size(x::SymbolicView) = size(x.v)
Base.firstindex(x::SymbolicView) = firstindex(x.v)
Base.lastindex(x::SymbolicView) = lastindex(x.v)
Base.size(x::SymbolicView) = size(x.v)
Base.iterate(x::SymbolicView) = iterate(x.v)
Base.iterate(x::SymbolicView, state) = iterate(x.v, state)
Base.length(x::SymbolicView{N}) = N
Base.IndexStyle(::Type{SymbolicView}) = IndexLinear()
Base.getindex(x::SymbolicView, index) = s.v[_sym_to_int(x, index)]
Base.setindex!(x::SymbolicView, value, index) = s.v[_sym_to_int(x, index)] = value
function _sym_to_int(x::SymbolicView, sym::Symbol)
    idx = findfirst(isequal(sym), x.syms)
    isnothing(idx) && throw(ArgumentError("SymbolError: try to access SymbolicView($(x.syms)) with symbol $sym"))
    idx
end
_sym_to_int(x::SymbolicView, idx::Int) == idx
_sym_to_int(x::SymbolicView, idx) == _sym_to_int.(Ref(x), idx)
