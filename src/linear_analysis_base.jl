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
@kwdef struct NetworkDescriptorSystem{MT,AT,BT,CT,DT,ST,IST,OST}
    M::MT = LinearAlgebra.I
    A::AT
    B::BT = nothing
    C::CT = nothing
    D::DT = nothing
    sym::ST = nothing
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
        new{typeof(M), typeof(A), typeof(B), typeof(C), typeof(D), typeof(sym), typeof(insym), typeof(outsym)}(M, A, B, C, D, sym, insym, outsym)
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
                constraints = maybe_plural(dim(sys) - sum(sys.M), "constraint")[2]
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
            if insym(sys) isa AbstractVector
                _insymstr = "[" * join(string.(insym(sys)), ", ") * "]"
            else
                _insymstr = string(insym(sys))
            end
            if outsym(sys) isa AbstractVector
                _outsymstr = "[" * join(string.(outsym(sys)), ", ") * "]"
            else
                _outsymstr = string(outsym(sys))
            end
            push!(strvec, "$(ni) &$(_inputs): &$(_insymstr)")
            push!(strvec, "$(no) &$(_outputs): &$(_outsymstr)")
            print_treelike(io, align_strings(strvec))
        end
    end
end

####
#### LTI composition operations
####
"""
    _blkdiag(A, B, nothing, C,...)
    _blkdiag((A, n1), (B, n2), (nothing, _), (C, n3),...)

Block-diagonal concatenation of matrices. Ignore nothing.
Also works with UniformScaling if the size is provided (second signature).
"""
function _blkdiag(mats::Union{AbstractMatrix, Nothing}...)
    _mats = unrolled_map(mats) do mat
        isnothing(mat) ? (nothing, (0, 0)) : (mat, size(mat))
    end
    _blkdiag(_mats...)
end
_blkdiag(mats::Tuple{Nothing, <:Any}...) = nothing
_blkdiag(mats::Nothing...) = nothing
function _blkdiag(mats::Tuple...)
    I = Int[]
    J = Int[]
    V = cachetype(unrolled_map(first, mats)...)[]

    nrows, ncols = _collect_ijv_recursive!(I, J, V, mats)

    sparse(I, J, V, nrows, ncols)
end
function _collect_ijv_recursive!(I, J, V, mats::Tuple, roff=0, coff=0)
    isempty(mats) && return (roff, coff)
    mat, (nrows, ncols) = first(mats)
    if !isnothing(mat)
        if mat isa LinearAlgebra.UniformScaling
            for i in 1:nrows  # nrows == ncols for UniformScaling
                push!(I, i + roff)
                push!(J, i + coff)
                push!(V, mat.λ)
            end
        elseif SparseArrays.issparse(mat)
            @assert size(mat) == (nrows, ncols)
            _I, _J, _V = SparseArrays.findnz(mat)
            _I .+= roff
            _J .+= coff
            append!(I, _I)
            append!(J, _J)
            append!(V, _V)
        else
            @assert size(mat) == (nrows, ncols)
            @inbounds for matidx in eachindex(IndexCartesian(), mat)
                push!(I, matidx[1] + roff)
                push!(J, matidx[2] + coff)
                push!(V, mat[matidx])
            end
        end
        roff += nrows
        coff += ncols
    end
    _collect_ijv_recursive!(I, J, V, Base.tail(mats), roff, coff)
end

_cat_syms(syms...) = vcat(syms...)
_cat_syms(tups::Tuple{Nothing, <:Any}...) = nothing
function _cat_syms(tups::Tuple...)
    syms = unrolled_map(tups) do tup
        s, n = tup
        isnothing(s) ? fill(nothing, n) : s
    end
    vcat(syms...)
end

"""
    append(sys1, sys2, ...) -> NetworkDescriptorSystem

Block-diagonal (direct sum) of descriptor systems — independent systems sharing no
states or signals (MATLAB: `append`):

    M·ẋ = blkdiag(A1,A2) x + blkdiag(B1,B2) u
      y = blkdiag(C1,C2) x + blkdiag(D1,D2) u

`insym` and `outsym` are the concatenation of those of the subsystems.
"""
function append(systems::NetworkDescriptorSystem...)
    Ms = unrolled_map(s -> (s.M, size(s.A)), systems)
    As = unrolled_map(s -> s.A, systems)
    Bs = unrolled_map(s -> s.B, systems)
    Cs = unrolled_map(s -> s.C, systems)
    Ds = unrolled_map(s -> s.D, systems)
    syms = unrolled_map(s -> (s.sym, size(s.A, 1)), systems)
    insyms = unrolled_map(s -> (s.insym, indim(s)), systems)
    outsyms = unrolled_map(s -> (s.outsym, outdim(s)), systems)

    M = _blkdiag(Ms...)
    A = _blkdiag(As...)
    B = _blkdiag(Bs...)
    C = _blkdiag(Cs...)
    D = _blkdiag(Ds...)
    sym = _cat_syms(syms...)
    insym = _cat_syms(insyms...)
    outsym = _cat_syms(outsyms...)

    NetworkDescriptorSystem(; M, A, B, C, D, sym, insym, outsym)
end

"""
    sys2 * sys1 -> NetworkDescriptorSystem

Series connection: output of `sys1` feeds input of `sys2` (MATLAB: `series`).
Requires `outdim(sys1) == indim(sys2)`.

The combined state is `[x1; x2]`:

    A = [A1        0 ]    B = [B1      ]
        [B2·C1    A2 ]        [B2·D1   ]
    C = [D2·C1    C2 ]    D = D2·D1
"""
function Base.:*(sys2::NetworkDescriptorSystem, sys1::NetworkDescriptorSystem)
    if outdim(sys1) != indim(sys2)
        throw(DimensionMismatch("outdim(sys1)=$(outdim(sys1)) must equal indim(sys2)=$(indim(sys2))"))
    end
    n1, n2 = dim(sys1), dim(sys2)
    M = _blkdiag((sys1.M, size(sys1.A)), (sys2.M, size(sys2.A)))
    A = [sys1.A              zeros(n1, n2);
         sys2.B * sys1.C     sys2.A       ]
    B = [sys1.B;
         sys2.B * sys1.D]
    C = [sys2.D * sys1.C    sys2.C]
    D = sys2.D * sys1.D
    sym = _cat_syms((sys1.sym, n1), (sys2.sym, n2))
    return NetworkDescriptorSystem(; M, A, B, C, D, sym,
        insym=sys1.insym, outsym=sys2.outsym)
end

"""
    k * sys -> NetworkDescriptorSystem
    sys * k -> NetworkDescriptorSystem

Scalar gain: scales the output of `sys` by `k`.

    C_new = k * C,  D_new = k * D

`A`, `B`, and `M` are unchanged. Unary negation `-sys` is `(-1) * sys`.
"""
function Base.:*(k::Number, sys::NetworkDescriptorSystem)
    NetworkDescriptorSystem(; M=sys.M, A=sys.A, B=sys.B, C=k*sys.C, D=k*sys.D,
        sym=sys.sym, insym=sys.insym, outsym=sys.outsym)
end
Base.:*(sys::NetworkDescriptorSystem, k::Number) = k * sys
Base.:-(sys::NetworkDescriptorSystem) = (-1) * sys
Base.:-(sys1::NetworkDescriptorSystem, sys2::NetworkDescriptorSystem) = sys1 + (-sys2)

"""
    Σ * sys -> NetworkDescriptorSystem

Left-multiply (output transformation): `C_new = Σ*C`, `D_new = Σ*D`.
`A`, `B`, `M` are unchanged. `outsym` is replaced by generic symbols.
"""
function Base.:*(Σ::AbstractMatrix, sys::NetworkDescriptorSystem)
    n_out = size(Σ, 1)
    NetworkDescriptorSystem(; M=sys.M, A=sys.A, B=sys.B, C=Σ*sys.C, D=Σ*sys.D,
        sym=sys.sym, insym=sys.insym,
        outsym=[Symbol("out_", i) for i in 1:n_out])
end

"""
    sys * Γ -> NetworkDescriptorSystem

Right-multiply (input transformation): `B_new = B*Γ`, `D_new = D*Γ`.
`A`, `C`, `M` are unchanged. `insym` is replaced by generic symbols.
"""
function Base.:*(sys::NetworkDescriptorSystem, Γ::AbstractMatrix)
    n_in = size(Γ, 2)
    NetworkDescriptorSystem(; M=sys.M, A=sys.A, B=sys.B*Γ, C=sys.C, D=sys.D*Γ,
        sym=sys.sym,
        insym=[Symbol("in_", i) for i in 1:n_in],
        outsym=sys.outsym)
end

"""
    sys1 + sys2 -> NetworkDescriptorSystem

Parallel connection: same inputs, outputs summed (MATLAB: `parallel`).
Requires matching input and output dimensions.

The combined state is `[x1; x2]`:

    A = blkdiag(A1,A2)    B = [B1; B2]
    C = [C1  C2]           D = D1 + D2
"""
function Base.:+(sys1::NetworkDescriptorSystem, sys2::NetworkDescriptorSystem)
    if indim(sys1) != indim(sys2)
        throw(DimensionMismatch("indim(sys1)=$(indim(sys1)) must equal indim(sys2)=$(indim(sys2))"))
    end
    if outdim(sys1) != outdim(sys2)
        throw(DimensionMismatch("outdim(sys1)=$(outdim(sys1)) must equal outdim(sys2)=$(outdim(sys2))"))
    end
    n1, n2 = dim(sys1), dim(sys2)
    M = _blkdiag((sys1.M, size(sys1.A)), (sys2.M, size(sys2.A)))
    A = _blkdiag(sys1.A, sys2.A)
    B = vcat(sys1.B, sys2.B)
    C = hcat(sys1.C, sys2.C)
    D = sys1.D + sys2.D
    sym = _cat_syms((sys1.sym, n1), (sys2.sym, n2))
    return NetworkDescriptorSystem(; M, A, B, C, D, sym,
        insym=sys1.insym, outsym=sys1.outsym)
end
