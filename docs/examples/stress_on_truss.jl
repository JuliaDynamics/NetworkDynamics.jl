#=
# Stress on Truss
In this exampe we'll simulate the time evolution of a truss structure consisting
of joints and stiff springs.

![truss animation](truss.mp4)

The mathematical model is quite simple: the vertices are point masses with
positions `x` and `y` and velocities `vx` and `vy`. The edges are relativly stiff springs.
```@math
\dot{x} &= vx
\dot{y} &= vy
\dot{vx} &= (\sum F_x - γ*vx) / M
\dot{vy} &= (\sum F_y - γ*vy) / M - g
```
The vertices cannot absorb any torque, so the beams only exert forces in the direction of the beam.

```@math
|F| = K\cdot(L - |Δd|)
```

We start by loading the necessary packages.
=#
using NetworkDynamics
using OrdinaryDiffEqTsit5
using Graphs
using GraphMakie
using LinearAlgebra: norm, ⋅
using Printf
using CairoMakie
CairoMakie.activate!()

#=
## Definition of the dynamical system

We need 3 models:
- a fixed vertex which cannot change its position,
- a free vertex which can move,
- a beam which connects two vertices.
=#

function fixed_g(pos, u, p, t)
    pos .= p
end
vertex_fix = VertexFunction(g=fixed_g, psym=[:xfix, :yfix], outsym=[:x, :y], ff=NoFeedForward())
#-

function free_f(dx, x, Fsum, (M, γ, g), t)
    v = view(x, 1:2)
    dx[1:2] .= (Fsum .- γ .* v) ./ M
    dx[2] -= g
    dx[3:4] .= v
    nothing
end
vertex_free = VertexFunction(f=free_f, g=3:4, sym=[:vx=>0, :vy=>0, :x, :y],
                             psym=[:M=>10, :γ=>200, :g=>9.81], insym=[:Fx, :Fy])
#=
For the edge, we want to dosomething special. Later on, we want to color the edges according to the force they exert.
Therefore, are interested in the absolut force rather than just the force vector.
NetworkDynamics allows you to  define so called `Observed` functions, which can recover additional states,
so called observed, after the simulations
=#

function edge_g!(F, pos_src, pos_dst, (K, L), t)
    dx = pos_dst[1] - pos_src[1]
    dy = pos_dst[2] - pos_src[2]
    d = sqrt(dx^2 + dy^2)
    Fabs = K * (L - d)
    F[1] = Fabs * dx / d
    F[2] = Fabs * dy / d
    nothing
end
function observedf(obsout, u, pos_src, pos_dst, (K, L), t)
    dx = pos_dst[1] .- pos_src[1]
    dy = pos_dst[2] .- pos_src[2]
    d = sqrt(dx^2 + dy^2)
    obsout[1] = K * (L - d)
    nothing
end
beam = EdgeFunction(g=AntiSymmetric(edge_g!), psym=[:K=>0.5e6, :L], outsym=[:Fx, :Fy], obsf=observedf, obssym=[:Fabs])

#=
Set up graph topology and initial positions.
=#
N = 5
dx = 1.0
shift = 0.2
g = SimpleGraph(2*N + 1)
for i in 1:N
    add_edge!(g, i, i+N); add_edge!(g, i, i+N)
    if i < N
        add_edge!(g, i+1, i+N); add_edge!(g, i, i+1); add_edge!(g, i+N, i+N+1)
    end
end
add_edge!(g, 2N, 2N+1)
pos0 = zeros(Point2f, 2N + 1)
pos0[1:N] = [Point((i-1)dx,0) for i in 1:N]
pos0[N+1:2*N] = [Point(i*dx + shift, 1) for i in 1:N]
pos0[2N+1] = Point(N*dx + 1, -1)
fixed = [1,4] # set fixed vertices
nothing #hide

#=
Now can collect the vertex models.
=#
verts = VertexFunction[vertex_free for i in 1:nv(g)]
for i in fixed
    verts[i] = vertex_fix # use the fixed vertex for the fixed points
end
nw = Network(g, verts, beam)
#=
In order to simulate the system we need to initialize the state and parameter vectors.

Some states and parameters are shared between all vertices/edges. Those have been allready set
in their constructors. The free symbols are
- `x` and `y` for the position of the free vertices,
- `xfix` and `yfix` for the position of the fixed vertices,
- `L` for the length of the beams.
=#
s = NWState(nw)
## set x/y and xfix/yfix
for i in eachindex(pos0, verts)
    if i in fixed
        s.p.v[i, :xfix] = pos0[i][1]
        s.p.v[i, :yfix] = pos0[i][2]
    else
        s.v[i, :x] = pos0[i][1]
        s.v[i, :y] = pos0[i][2]
    end
end
## set L for edges
for (i,e) in enumerate(edges(g))
    s.p.e[i, :L] = norm(pos0[src(e)] - pos0[dst(e)])
end
nothing #hide

#=
Lastly there is a special vertex at the end of the truss which has a higher mass.
=#
s.p.v[11, :M] = 200
s.p.v[11, :γ] = 100
nothing #hide

#=
No we have everything ready to build the ODEProblem and simulate the system.
=#
tspan = (0.0, 12.0)
prob = ODEProblem(nw, uflat(s), tspan, pflat(s))
sol  = solve(prob, Tsit5())
nothing #hide

#=
## Plot the solution

To plot the solution we want to make use of the `GraphMakie` package.
=#
fig = Figure(size=(1000,550));
fig[1,1] = title = Label(fig, "Stress on truss", fontsize=30)
title.tellwidth = false

fig[2,1] = ax = Axis(fig)

## get the maximum force during the simulation to get the color scale
(fmin, fmax) = 0.3 .*extrema(Iterators.flatten(sol(sol.t, idxs=eidxs(nw, :, :Fabs))))
nlabels = ["" for i in 1:nv(g)]
nlabels[fixed] .= "Δ"
elabels = ["edge $i" for i in 1:ne(g)]
p = graphplot!(ax, g;
               edge_width=4.0,
               node_size=3*sqrt.(try s.p.v[i, :M] catch; 10.0 end for i in 1:nv(g)),
               nlabels=nlabels,
               nlabels_align=(:center,:top),
               nlabels_fontsize=30,
               elabels=elabels,
               elabels_side=Dict(ne(g) => :right),
               edge_color=[0.0 for i in 1:ne(g)],
               edge_attr=(colorrange=(fmin,fmax),
                          colormap=:diverging_bkr_55_10_c35_n256))
hidespines!(ax); hidedecorations!(ax);
limits!(ax, -0.1, pos0[end][1]+0.3, pos0[end][2]-0.5, 1.15)

ax.aspect = DataAspect();

## draw colorbar
fig[3,1] = cb = Colorbar(fig, p.plots[1].plots[1], label = "Axial force", vertical=false)

T = tspan[2]
fps = 30
trange = range(0.0, sol.t[end], length=Int(T * fps))
record(fig, "truss.mp4", trange; framerate=fps) do t
    title.text = @sprintf "Stress on truss (t = %.2f )" t
    s_at_t = NWState(sol, t)
    for i in eachindex(pos0)
        p[:node_pos][][i] = (s_at_t.v[i, :x], s_at_t.v[i, :y])
    end
    notify(p[:node_pos])
    load = s_at_t.e[:, :Fabs]
    p.edge_color[] = load
    p.elabels = [@sprintf("%.0f", l) for l in load]
    fig
end
