abstract type ComponentCallback{C,A} end

struct ContinousComponentCallback{C,A} <: ComponentCallback{C,A}
    condition::C
    affect::A
    kwargs::KW
end

struct VectorContinousComponentCallback{C,A} <: ComponentCallback{C,A}
    condition::C
    affect::A
    len::Int
    kwargs::KW
end

struct DiscreteComponentCallback{C,A} <: ComponentCallback{C,A}
    condition::C
    affect::A
    kwargs::KW
end

struct CallbackBatch{T,C,A,DIM,PDIM}
    idxs::Vector{Int}
    urs::Vector{UnitRange{Int}}
    usyms::NTuple{DIM,Symbol}
    prs::Vector{UnitRange{Int}}
    psyms::NTuple{PDIM,Symbol}
    condition::C
    affect::A
    "Map outidx -> cb idx"
    reversemap::Vector{Int}
end

function symbolic_view_u(cbbatch::CallbackBatch, u, i)
    _u = view(u, cbbatch.urs[i])
    SymbolicView(_u, cbbatch.usyms)
end
function symbolic_view_p(cbbatch::CallbackBatch, p, i)
    _p = view(p, cbbatch.prs[i])
    SymbolicView(_p, cbbatch.psyms)
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

function batch_condition(cbbatch)
    (out, u, t, integrator) -> begin
        for i in 1:length(cbbatch)
            _u = symbolic_view_u(cbbatch, u, i)
            _p = symbolic_view_p(cbbatch, p, i)

            if cbtype(cbbatch) <: VectorContinousComponentCallback
                _out = batch.condition(_u, _p, t)
            elseif cbtype(cbbatch) <: ContinousComponentCallback
                @views _out = out[outrange(cbbatch, i)]
                cbbatch.condition(_out, _u, _p, t)
            end
        end
        nothing
    end
end
function batch_affect(cbbatch)
    (integrator, outidx) -> begin
        i = cbbatch.reversmap[outidx]
        _u = symbolic_view_u(cbbatch, u, i)
        _p = symbolic_view_p(cbbatch, p, i)

        uhash = hash(_u.v)
        phash = hash(_p.v)
        cbbatch.affect(_u, _p, integrator.t)
        pchanged = hash(_p.v) != phash
        uchanged = hash(_u.v) != uhash

        (pchanged || uchanged) && auto_dt_reset!(integrator)
        pchanged && save_parameters!(integrator)
    end
end
