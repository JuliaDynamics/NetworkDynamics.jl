# A tiny power grid in plain MTK: a powerflow network is solved with the sparse fixpoint solver
# and hands its interface values to the componentwise init of the dynamic network. Downstream
# packages do the same on every first run, so most of the generic code lands in the pkgimage.
using ModelingToolkitBase: @variables, @parameters, @named, System, t_nounits as t, D_nounits as Dt
using NetworkDynamics: Network, VertexModel, EdgeModel, NWState, VIndex, SymbolicIndex,
                       InitConstraint, InitFormula, set_jac_prototype!, find_fixpoint,
                       initialize_componentwise, interface_values, implicit_output

# voltage out, summed line currents in
function _wl_busbar(; name)
    ps = @parameters Vbase [initf_weak = 1.0]
    vars = @variables begin
        u_r(t), [guess=1, output=true]
        u_i(t), [guess=0, output=true]
        i_r(t), [input=true]
        i_i(t), [input=true]
        P(t)
        Q(t)
        u_kV(t)
    end
    eqs = [P ~ -u_r*i_r - u_i*i_i, Q ~ -u_i*i_r + u_r*i_i, u_kV ~ sqrt(u_r^2 + u_i^2)*Vbase]
    System(eqs, t, vars, ps; name)
end

# terminal variables of an injector, `i` is the current into the bus
_wl_terminal() = @variables u_r(t) u_i(t) i_r(t) [guess=0] i_i(t) [guess=0]

function _wl_swing(; name)
    term = _wl_terminal()
    u_r, u_i, i_r, i_i = term
    ps = @parameters M=0.005 D=0.1 V [guess=1, initf=sqrt(u_r^2 + u_i^2)] Pm [guess=1]
    vars = @variables θ(t) [guess=0, initf=atan(u_i, u_r)] ω(t)=0 Pel(t)
    eqs = [Dt(θ) ~ ω,
           Dt(ω) ~ (Pm - D*ω - Pel)/M,
           Pel ~ u_r*i_r + u_i*i_i,
           u_r ~ V*cos(θ),
           u_i ~ V*sin(θ)]
    System(eqs, t, vcat(term, vars), ps; name)
end

function _wl_pqload(; name)
    term = _wl_terminal()
    u_r, u_i, i_r, i_i = term
    ps = @parameters P=-0.5 Q=-0.1
    System([i_r ~ (P*u_r + Q*u_i)/(u_r^2 + u_i^2),
            i_i ~ (P*u_i - Q*u_r)/(u_r^2 + u_i^2)], t, term, ps; name)
end

function _wl_slack(; name)
    term = _wl_terminal()
    u_r, u_i, i_r, i_i = term
    ps = @parameters V=1 δ=0
    System([u_r ~ V*cos(δ), u_i ~ V*sin(δ)], t, term, ps; name)
end

function _wl_bus(injector; name, vidx)
    @named busbar = _wl_busbar()
    eqs = [injector.u_r ~ busbar.u_r, injector.u_i ~ busbar.u_i,
           0 ~ busbar.i_r + injector.i_r, 0 ~ busbar.i_i + injector.i_i]
    sys = System(eqs, t; name, systems=[busbar, injector])
    VertexModel(sys, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; vidx)
end

# a bus without injector, the currents only constrain the voltage implicitly
function _wl_junction(; name, vidx)
    @named busbar = _wl_busbar()
    sys = System([busbar.i_r ~ implicit_output(busbar.u_r),
                  busbar.i_i ~ implicit_output(busbar.u_i)], t; name, systems=[busbar])
    VertexModel(sys, [:busbar₊i_r, :busbar₊i_i], [:busbar₊u_r, :busbar₊u_i]; vidx)
end

function _wl_lineend(; name)
    @variables u_r(t) [input=true] u_i(t) [input=true] i_r(t) [output=true] i_i(t) [output=true] P(t)
    System([P ~ u_r*i_r + u_i*i_i], t; name)
end

function _wl_piline(; name, src, dst)
    @named src_end = _wl_lineend()
    @named dst_end = _wl_lineend()
    @parameters R=0.01 X=0.1 B=0.02
    @variables i_r(t) i_i(t)
    Δr, Δi = src_end.u_r - dst_end.u_r, src_end.u_i - dst_end.u_i
    eqs = [i_r ~ (R*Δr + X*Δi)/(R^2 + X^2),
           i_i ~ (R*Δi - X*Δr)/(R^2 + X^2),
           src_end.i_r ~ -i_r + B/2*src_end.u_i,
           src_end.i_i ~ -i_i - B/2*src_end.u_r,
           dst_end.i_r ~ i_r + B/2*dst_end.u_i,
           dst_end.i_i ~ i_i - B/2*dst_end.u_r]
    sys = System(eqs, t; name, systems=[src_end, dst_end])
    EdgeModel(sys, [:src_end₊u_r, :src_end₊u_i], [:dst_end₊u_r, :dst_end₊u_i],
                   [:src_end₊i_r, :src_end₊i_i], [:dst_end₊i_r, :dst_end₊i_i]; src, dst)
end

# the parallel line gives a `ComponentGraph`, as in most real grids
_wl_lines() = [_wl_piline(; name=:l12a, src=1, dst=2),
               _wl_piline(; name=:l12b, src=1, dst=2),
               _wl_piline(; name=:l13, src=1, dst=3),
               _wl_piline(; name=:l23, src=2, dst=3)]

function _wl_run()
    pfnw = Network([_wl_bus(_wl_slack(; name=:slack); name=:bus1, vidx=1),
                    _wl_junction(; name=:bus2, vidx=2),
                    _wl_bus(_wl_pqload(; name=:load); name=:bus3, vidx=3)], _wl_lines())
    set_jac_prototype!(pfnw)
    pfs = find_fixpoint(pfnw, NWState(pfnw; guess=true); alg=nothing)

    nw = Network([_wl_bus(_wl_swing(; name=:swing); name=:bus1, vidx=1),
                  _wl_junction(; name=:bus2, vidx=2),
                  _wl_bus(_wl_pqload(; name=:load); name=:bus3, vidx=3)], _wl_lines())
    # kwargs in the order PowerDynamics passes them, the kwcall is specialized on it
    s0 = initialize_componentwise(nw; default_overrides=interface_values(pfs),
                                  additional_initconstraint=Dict{SymbolicIndex,Vector{InitConstraint}}(),
                                  additional_initformula=Dict{SymbolicIndex,Vector{InitFormula}}(),
                                  verbose=false, subverbose=false)
    s0[VIndex(1, :busbar₊P)]
end
_wl_run()
