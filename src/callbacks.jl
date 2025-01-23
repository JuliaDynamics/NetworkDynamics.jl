abstract type ComponentCallback{C,A} end

struct ComponentCondition{C,DIM,PDIM}
    f::C
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
    function ComponentCondition(f, sym, psym)
        new{typeof(f), length(sym), length(psym)}(f, Tuple(sym), Tuple(psym))
    end
end

struct ComponentAffect{A,DIM,PDIM}
    f::A
    sym::NTuple{DIM,Symbol}
    psym::NTuple{PDIM,Symbol}
    function ComponentAffect(f, sym, psym)
        new{typeof(f), length(sym), length(psym)}(f, Tuple(sym), Tuple(psym))
    end
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
    kwargs::NamedTuple
end
function ContinousComponentCallback(condition, affect; kwargs...)
    ContinousComponentCallback(condition, affect, NamedTuple(kwargs))
end

struct VectorContinousComponentCallback{C,A,CDIM,CPDIM,ADIM,APDIM} <: ComponentCallback{C,A}
    condition::ComponentCondition{C,CDIM,CPDIM}
    affect::ComponentAffect{A,ADIM,APDIM}
    len::Int
    kwargs::NamedTuple
end
function VectorContinousComponentCallback(condition, affect, len; kwargs...)
    VectorContinousComponentCallback(condition, affect, len, NamedTuple(kwargs))
end

struct CallbackBatch{T<:ComponentCallback,C,A,ST<:SymbolicIndex}
    nw::Network
    components::Vector{ST}
    callbacks::Vector{T}
    sublen::Int # length of each callback
    condition::C
    affect::A
end
function CallbackBatch(nw, components, callbacks)
    if !isconcretetype(eltype(components))
        components = Vector{typeof(first(components))}(components)
    end
    if !isconcretetype(eltype(callbacks))
        callbacks = Vector{typeof(first(callbacks))}(callbacks)
    end
    sublen = eltype(callbacks) <: ContinousComponentCallback ? 1 : first(callbacks).len
    condition = first(callbacks).condition.f
    affect = first(callbacks).affect.f
    CallbackBatch(nw, components, callbacks, sublen, condition, affect)
end

Base.length(cbb::CallbackBatch) = length(cbb.callbacks)

cbtype(cbb::CallbackBatch{T}) where {T} = T

condition_dim(cbb)  = first(cbb.callbacks).condition.sym  |> length
condition_pdim(cbb) = first(cbb.callbacks).condition.psym |> length
affect_dim(cbb)  = first(cbb.callbacks).affect.sym  |> length
affect_pdim(cbb) = first(cbb.callbacks).affect.psym |> length

condition_urange(cbb, i) = (1 + (i-1)*condition_dim(cbb))  : i*condition_dim(cbb)
condition_prange(cbb, i) = (1 + (i-1)*condition_pdim(cbb)) : i*condition_pdim(cbb)
affect_urange(cbb, i) = (1 + (i-1)*affect_dim(cbb) ) : i*affect_dim(cbb)
affect_prange(cbb, i) = (1 + (i-1)*affect_pdim(cbb)) : i*affect_pdim(cbb)

condition_outrange(cbb, i) = (1 + (i-1)*cbb.sublen) : i*cbb.sublen

cbidx_from_outidx(cbb, outidx) = div(cbb.sublen, outidx) + 1

function collect_c_or_a_indices(cbb::CallbackBatch, c_or_a, u_or_p)
    sidxs = SymbolicIndex[]
    for (component, cb) in zip(cbb.components, cbb.callbacks)
        syms = getproperty(getproperty(cb, c_or_a), u_or_p)
        symidxtype = if component isa VIndex
            u_or_p == :sym ? VIndex : VPIndex
        else
            u_or_p ==:sym ? EIndex : EPIndex
        end
        sidx = collect(symidxtype(component.compidx, syms))
        append!(sidxs, sidx)
    end
    sidxs
end

function collect_callbackbatches(nw)
    components = SymbolicIndex[]
    callbacks = ComponentCallback[]
    for (i, v) in pairs(nw.im.vertexm)
        has_callback(v) || continue
        for cb in get_callbacks(v)
            push!(components, VIndex(i, nothing))
            push!(callbacks, cb)
        end
    end
    for (i, v) in pairs(nw.im.edgem)
        has_callback(v) || continue
        for cb in get_callbacks(v)
            push!(components, EIndex(i, nothing))
            push!(callbacks, cb)
        end
    end

    idx_per_type = _find_identical(callbacks, 1:length(components))
    batches = CallbackBatch[]
    for typeidx in idx_per_type
        batchcomps = components[typeidx]
        batchcbs = callbacks[typeidx]
        cbs = CallbackBatch(nw, batchcomps, batchcbs)
        push!(batches, cbs)
    end
    return batches
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

function batch_condition(cbb)
    usymidxs = collect_c_or_a_indices(cbb, :condition, :sym)
    psymidxs = collect_c_or_a_indices(cbb, :condition, :psym)
    ucache = DiffCache(zeros(length(usymidxs)), 12)
    obsf = SII.observed(cbb.nw, usymidxs)
    pidxs = SII.parameter_index.(Ref(cbb.nw), psymidxs)

    (out, u, t, integrator) -> begin
        us = PreallocationTools.get_tmp(ucache, u)
        obsf(u, integrator.p, t, us)

        for i in 1:length(cbb)
            # symbolic view into u
            uv = view(us, condition_urange(cbb, i))
            _u = SymbolicView(uv, cbb.callbacks[i].condition.sym)

            # symbolic view into p
            pidxsv = view(pidxs, condition_prange(cbb, i))
            pv = view(integrator.p, pidxsv)
            _p = SymbolicView(pv, cbb.callbacks[i].condition.psym)

            if cbtype(cbb) <: ContinousComponentCallback
                oidx = only(condition_outrange(cbb, i))
                out[oidx] = cbb.condition(_u, _p, t)
            elseif cbtype(cbb) <: VectorContinousComponentCallback
                @views _out = out[condition_outrange(cbb, i)]
                cbb.condition(_out, _u, _p, t)
            end
        end
        nothing
    end
end
function batch_affect(cbb)
    usymidxs = collect_c_or_a_indices(cbb, :affect, :sym)
    psymidxs = collect_c_or_a_indices(cbb, :affect, :psym)
    uidxs = SII.variable_index.(Ref(cbb.nw), usymidxs)
    pidxs = SII.parameter_index.(Ref(cbb.nw), psymidxs)

    (integrator, outidx) -> begin
        uidxsv = view(uidxs, affect_urange(cbb, i))
        uv = view(integrator.u, uidxsv)
        _u = SymbolicView(uv, cbb.callbacks[i].affect.sym)

        pidxsv = view(pidxs, affect_prange(cbb, i))
        pv = view(integrator.p, pidxsv)
        _p = SymbolicView(pv, cbb.callbacks[i].affect.psym)

        i = cbidx_from_outidx(cbb, outidx)

        uhash = hash(uv)
        phash = hash(pv)
        cbb.affect(_u, _p, integrator.t, get_ctx(cbb, i))
        pchanged = hash(pv) != phash
        uchanged = hash(uv) != uhash

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
Base.iterate(x::SymbolicView) = iterate(x.v)
Base.iterate(x::SymbolicView, state) = iterate(x.v, state)
Base.length(x::SymbolicView{N}) where {N} = N
Base.IndexStyle(::Type{SymbolicView}) = IndexLinear()
Base.getindex(x::SymbolicView, index) = x.v[_sym_to_int(x, index)]
Base.setindex!(x::SymbolicView, value, index) = x.v[_sym_to_int(x, index)] = value
function _sym_to_int(x::SymbolicView, sym::Symbol)
    idx = findfirst(isequal(sym), x.syms)
    isnothing(idx) && throw(ArgumentError("SymbolError: try to access SymbolicView($(x.syms)) with symbol $sym"))
    idx
end
_sym_to_int(x::SymbolicView, idx::Int) = idx
_sym_to_int(x::SymbolicView, idx) = _sym_to_int.(Ref(x), idx)
