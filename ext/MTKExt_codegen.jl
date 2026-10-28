# Builds f or g. For `Dt(x) ~ -y*k`, `y ~ sin(x)` and `o ~ ifelse(x > 0, y, 0)` the generated g
# looks roughly like this:
#
#     function g(out, u)
#         x = u[1]                    # argument elements become named locals
#         ##cse#1 = 0 < x
#         o = if ##cse#1
#             ##cse#2 = sin(x)        # only computed in the branch which needs it
#             y = ##cse#2
#             y
#         else
#             0
#         end
#         out[1] = o                  # all output writes at the end
#         out
#     end
function build_fg(eqs, obs_subs, args; expression)
    assignments, outs = _get_formulas(eqs, obs_subs)
    _build_bound(assignments, outs, args; cse=true, expression)
end

# The obsf takes an extra mask argument which says which observed assignments to compute, see
# `DependencyAwareObsf`. No CSE here: it would move work shared by two observables out of their
# mask guards. With the observables `a ~ 2x` and `b ~ a + k`:
#
#     function obsf(out, u, in, p, t, mask)
#         a = if mask[1]              # every observed assignment has its own bit
#             x = u[1]
#             2x
#         else
#             NaN
#         end
#         b = if mask[2]              # set together with mask[1], as b needs a
#             k = p[1]
#             a + k
#         else
#             NaN
#         end
#         out[1] = a
#         out[2] = b
#         out
#     end
function build_obsf(obsstates, obs_subs, allinputs, args; expression)
    obsmask = unwrap(Symbolics.variable(:ˍ₋obsmask))
    assignments, columns, needs_input = _get_masked_obsformulas(obsstates, obs_subs, allinputs, obsmask)
    f_ip = _build_bound(assignments, obsstates, (args..., obsmask); cse=false, expression)
    NetworkDynamics.DependencyAwareObsf(f_ip, columns, needs_input)
end

# The assignments are a single `Let` around the output writes, and every argument element is bound
# to a named local (`state1 = ˍ₋arg1[1]`). `sink_conv` then moves each assignment down to where it
# is first needed.
function _build_bound(assignments, outs, args; cse, expression)
    wrap = fex -> Func(map(_with_bindings, fex.args), fex.kwargs, Let(assignments, fex.body, false), fex.pre)
    _, f_ip = build_function(outs, args...; cse, wrap_code=(wrap, wrap), conv=sink_conv, expression)
    f_ip
end
_with_bindings(arg) = arg
_with_bindings(d::DestructuredArgs) = DestructuredArgs(d.elems, d.inds, d.name, d.inbounds, true)

function _get_formulas(eqs, obs_subs)
    obsdeps = _collect_deps_on_obs([eq.rhs for eq in eqs], obs_subs)
    obs_assignments = [Assignment(k, v) for (k,v) in obs_subs if k ∈ obsdeps]

    # implicit equations do not become assignments (lhs is 0), so we use the rhs directly for output.
    # An output which another output reads is an observable too, it keeps its observable assignment.
    obskeys = Set(a.lhs for a in obs_assignments)
    eqs_assignments = [Assignment(eq.lhs, eq.rhs) for eq in eqs
                          if !isequal(eq.lhs, eq.rhs) && !isequal(unwrap_const(eq.lhs), 0) && eq.lhs ∉ obskeys]
    # implicit eqs have 0 lhs, explicit eqs are already captured via assignment, output their lhs
    outs = [isequal(unwrap_const(eq.lhs), 0) ? eq.rhs : eq.lhs for eq in eqs]

    # integer powers become plain products, so no `pow` call is left in the generated code
    topow = SymbolicUtils.Postwalk(_power_to_product)
    assignments = [Assignment(a.lhs, topow(a.rhs)) for a in vcat(obs_assignments, eqs_assignments)]
    assignments, map(topow, outs)
end

# Like `_get_formulas` for the observables, but every assignment is guarded by `mask`.
# Also returns what `DependencyAwareObsf` needs: per observable, the assignments it needs and
# whether one of them reads an input.
function _get_masked_obsformulas(obsstates, obs_subs, allinputs, mask)
    # filter out obs_subs which are not needed for any observable
    deps = _collect_deps_on_obs(obsstates, obs_subs)
    assigned = [k for k in keys(obs_subs) if k ∈ deps]
    topow = SymbolicUtils.Postwalk(_power_to_product)
    assignments = map(enumerate(assigned)) do (j, k)
        Assignment(k, ifelse(_obs_required(mask, j), topow(obs_subs[k]), NaN))
    end

    # position of each assigned symbol, also its bit in the mask
    assignment_index = Dict(k => j for (j, k) in enumerate(assigned))
    assignment_rhs_syms = [_all_rhs_symbols(obs_subs[k]) for k in assigned]

    # for each assignment, mark itself and every assignment it needs. The assignments are
    # topologically sorted, so the columns of its dependencies are already complete.
    assignment_columns = Vector{BitVector}(undef, length(assigned))
    for (j, syms) in enumerate(assignment_rhs_syms)
        col = falses(length(assigned))
        col[j] = true
        for s in syms
            haskey(assignment_index, s) || continue
            dep = assignment_index[s]
            dep < j || error("Observed assignments are not topologically sorted.")
            NetworkDynamics._or!(col, assignment_columns[dep])
        end
        assignment_columns[j] = col
    end
    obs_columns = [assignment_columns[assignment_index[s]] for s in obsstates]

    # an observable needs the inputs if any of its assignments reads one
    inputs = Set{ST}(allinputs) # including ext inputs
    assignment_reads_input = BitVector(any(∈(inputs), syms) for syms in assignment_rhs_syms)
    obs_needs_input = BitVector(any(assignment_reads_input .& col) for col in obs_columns)

    assignments, obs_columns, obs_needs_input
end
_obs_required(::Nothing, j) = true
_obs_required(mask, j) = @inbounds mask[j]
Symbolics.@register_symbolic _obs_required(mask, j::Int)::Bool

function _power_to_product(ex)
    iscall(ex) && operation(ex) === (^) || return ex
    x, n = arguments(ex)
    n = unwrap_const(n)
    n isa Integer && n >= 2 || return ex
    _power_product(x, n, SymbolicUtils.symtype(ex))
end
# x^n by repeated squaring, so x^100 takes 8 multiplications. A plain `x^n` would not do: CSE
# binds the exponent to a variable, and a power with a variable exponent is a slow generic call.
# A raw `Term` keeps `x*x` from collapsing back into `x^2`.
function _power_product(x, n, T)
    n == 1 && return x
    half = _power_product(x, n ÷ 2, T)
    sq = SymbolicUtils.Term{SymbolicUtils.vartype(x)}(*, [half, half]; type=T)
    isodd(n) ? SymbolicUtils.Term{SymbolicUtils.vartype(x)}(*, [sq, x]; type=T) : sq
end

# After CSE the generated code is a list of assignments followed by the output writes. Every
# assignment is pure and assigned once, so it can move down to just before its first user, and
# into an `if` arm when nothing outside that arm uses it. Unused assignments are dropped. The
# output writes stay together at the end: `out` might alias an input, so LLVM cannot move loads
# across a write.
sink_conv(ex, st) = _sink_assignments(Symbolics.toexpr(ex, st))
function _sink_assignments(fex)
    stmts = _flatten!(Any[], fex.args[2])
    k = something(findfirst(!_isassignment, stmts), length(stmts) + 1)
    defs, writes = stmts[1:k-1], stmts[k:end]
    # a name assigned twice can't be moved around, leave such code as it is
    allunique(a.args[1] for a in defs) || return fex
    # bottom-up, so the users of an assignment are already in place
    placed = Any[]
    for a in reverse(defs)
        s = a.args[1]
        if any(_mentions(s), writes)
            # needed at the top level, so it may not move into an `if` arm
            insert!(placed, something(findfirst(_mentions(s), placed), length(placed) + 1), a)
        elseif any(_mentions(s), placed)
            _place!(placed, a)
        end
    end
    Expr(:function, fex.args[1], Expr(:block, placed..., writes...))
end
_isassignment(ex) = Meta.isexpr(ex, :(=)) && ex.args[1] isa Symbol
_mentions(s) = ex -> _mentions(ex, s)
_mentions(ex, s) = ex === s
_mentions(ex::Expr, s) = any(_mentions(s), ex.args)
# The body is nested `begin` and `@inbounds begin` blocks, none of which opens a scope.
# `build_function` wraps the whole function in `@inbounds` again afterwards.
function _flatten!(stmts, ex)
    if Meta.isexpr(ex, :block)
        foreach(x -> _flatten!(stmts, x), ex.args)
    elseif Meta.isexpr(ex, :macrocall) && ex.args[1] === Symbol("@inbounds") && Meta.isexpr(ex.args[end], :block)
        _flatten!(stmts, ex.args[end])
    elseif !(ex isa LineNumberNode)
        push!(stmts, ex)
    end
    stmts
end
# Puts `a` right before its first user, or into that user's `if` arm if it is the only user.
function _place!(stmts, a)
    s = a.args[1]
    j = findfirst(_mentions(s), stmts)
    onlyuser = !any(_mentions(s), @view stmts[j+1:end])
    onlyuser && _sink_into!(stmts[j], a) && return
    insert!(stmts, j, a)
end
# Follows the single child which mentions the name down to an `if`. If the condition does not
# use it and only one arm does, the assignment goes into that arm. Loops, closures and the like
# are not entered, an assignment moved in there could run more than once.
function _sink_into!(ex, a)
    Meta.isexpr(ex, (:if, :block, :call, :(=))) || return false
    s = a.args[1]
    users = findall(_mentions(s), ex.args)
    length(users) == 1 || return false
    k = only(users)
    if ex.head === :if
        k == 1 && return false
        arm = Meta.isexpr(ex.args[k], :block) ? ex.args[k].args : Any[ex.args[k]]
        _place!(arm, a)
        ex.args[k] = Expr(:block, arm...)
        return true
    end
    _sink_into!(ex.args[k], a)
end
