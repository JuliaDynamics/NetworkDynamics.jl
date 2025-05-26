using NetworkDynamics
using Graphs
using Test

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

g = SimpleGraph([0 1 1 0 1;
                 1 0 1 1 0;
                 1 1 0 1 0;
                 0 1 1 0 1;
                 1 0 0 1 0])

vs = [Lib.swing_mtk() for _ in 1:5];
set_default!(vs[1], :Pmech, -1)
set_default!(vs[2], :Pmech, 1.5)
set_default!(vs[3], :Pmech, -1)
set_default!(vs[4], :Pmech, -1)
set_default!(vs[5], :Pmech, 1.5)

ls = [Lib.line_mtk() for _ in 1:7];
nw = Network(g, vs, ls)
s = NWState(nw)
p = NWParameter(nw)

repr(nw)
repr("text/plain", nw)

repr("text/plain", vs[1])
repr("text/plain", ls[1])


repr("text/plain", s)
repr("text/plain", p)
repr("text/plain", s.v)
repr("text/plain", s.e)
repr("text/plain", p.v)
repr("text/plain", p.e)

affect = ComponentAffect([], [:limit, :K]) do u, p, ctx
    nothing
end
condition = ComponentCondition([:P, :₋P, :srcθ], [:limit, :K]) do u, p, t
    nothing
end
cb1 = PresetTimeComponentCallback(1.0, affect)
repr("text/plain", cb1)
cb2 = ContinousComponentCallback(condition, affect)
repr("text/plain", cb2)
cb3 = DiscreteComponentCallback(condition, affect)
repr("text/plain", cb3)
set_callback!(ls[1], (cb1, cb2, cb3))
repr("text/plain", ls[1])
repr("text/plain", nw)
