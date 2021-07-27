using NDPrototype
using LightGraphs
using Random

@inline function diffusion_edge!(e,v_s,v_d,p,t)
    @inbounds for i in 1:length(e)
        e[i] = v_s[i] - v_d[i]
    end
    nothing
end

@inline function diffusion_vertex!(dv, v, agg, p, t)
    @inbounds for i in 1:length(dv)
        dv[i] = agg[i]
    end
    nothing
end

N = 10
g = SimpleGraph(watts_strogatz(N, 4, 0.5, seed=1))

# g = complete_graph(2)
vdim = 1
edim = 1
odevertex = ODEVertex(f = diffusion_vertex!, dim = vdim, pdim = 0)
staticedge = StaticEdge(f = diffusion_edge!, dim = edim, pdim = 0, coupling=AntiSymmetric())

nw = Network(g, odevertex, staticedge, verbose=true);
u = rand(dim(nw))
du = zeros(size(u)...)
p = Float64[]

@btime $nw($du, $u, $p, 0.0)

#############

odevertexv = [odevertex for i in 1:nv(g)];
staticedgev = [staticedge for i in 1:ne(g)];

nw = Network(g, odevertexv, staticedgev);

u = rand(dim(nw))
# u = [0.0, 1.0]
du = zeros(size(u)...)
p = Float64[]

nw(du, u, Float64[], 0.0)

#############

odevertex = [ODEVertex(f = (du,u,a,p,t)->diffusion_vertex!(du,u,a,p,i), dim = vdim, pdim = 0) for i in 1:nv(g)]
staticedge = [StaticEdge(f = (e,vs,vd,p,t)->diffusion_edge!(e,vs,vd,p,i), dim = edim, pdim = 0, coupling=AntiSymmetric()) for i in 1:ne(g)]

nw = Network(g, odevertex, staticedge);

u = rand(dim(nw))
# u = [0.0, 1.0]
du = zeros(size(u)...)
p = Float64[]

@btime $nw($du, $u, $p, 0.0)
cb = nw.nl.colorbatches[1]

#############
include("../benchmark/benchmark_utils.jl")

# diffusion test
N = 10
vertex = diffusion_vertex()
edge = diffusion_edge()
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end
# parallel
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    nd = Network(g, $vertex, $edge, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

# odeedge
@benchmark begin
    nd(dx, x0, nothing, 0.0)
end setup = begin
    g = watts_strogatz($N, 3, 0.8, seed=1)
    edge = diffusion_dedge()
    nd = Network(g, $vertex, edge, accdim=1, parallel=true)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, nothing, 0.0)
end

N = 10
g = watts_strogatz(N, 3, 0.8, seed=1)
edge = diffusion_dedge()
vertex = diffusion_vertex()
nd1 = Network(g, vertex, edge, accdim=1, parallel=true);
nd2 = Network(g, vertex, edge, accdim=1, parallel=false);
x0 = randn(dim(nd1));
dx1 = zero(x0);
dx2 = zero(x0);
nd1(dx1, x0, nothing, 0.0)
nd2(dx2, x0, nothing, 0.0)
@test dx1 ≈ dx2

# inhomgeneous network
function heterogeneous(N)
    rng = MersenneTwister(1)
    g = watts_strogatz(N, 3, 0.8, seed=1)
    edge = static_kuramoto_edge()
    vertex = [kuramoto_vertex_1d(), kuramoto_vertex_2d()]
    vertices = vertex[shuffle(rng, vcat([1 for _ in 1:N÷2], [2 for _ in 1:N÷2]))]
    p = vcat(randn(rng, nv(g)), randn(rng, ne(g)))
    (p, vertices, edge, g)
end

N = 10
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end
@benchmark begin
    nd(dx, x0, p, 0.0)
end setup = begin
    (p, v, e, g) = heterogeneous($N)
    nd = Network(g, v, e, parallel=true)
    @assert length(p) == pdim(nd)
    x0 = randn(dim(nd))
    dx = similar(x0)
    nd(dx, x0, p, 0.0) # call to init caches, we don't want to benchmark this
end

begin
    N = 10000
    (p, v, e, g) = heterogeneous(N)
    nd1 = Network(g, v, e, parallel=false)
    (p, v, e, g) = heterogeneous(N)
    nd2 = Network(g, v, e, parallel=true)
    x0 = randn(dim(nd1))
    dx1 = zeros(dim(nd1))
    dx2 = zeros(dim(nd1))
    @time nd1(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd1(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2

    @time nd2(dx1, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @time nd2(dx2, x0, p, 0.0) # call to init caches, we don't want to benchmark this
    @test dx1 ≈ dx2
    extrema(dx1 - dx2)

    @benchmark $nd1($dx1, $x0, $p, 0.0) # call to init caches, we don't want to benchmark this
    @benchmark $nd2($dx2, $x0, $p, 0.0) # call to init caches, we don't want to benchmark this
end

async_foo(i) = Threads.@spawn foo(i)

function foo(i)
    # @Threads.threads
    for j in 'a':'c'
        bar(i, j)
    end
end

function bar(i, j)
    sleep(.2)
    println("Thread=", Threads.threadid(), " i=$i j=$j")
end

function entry()
    async_foo(1)
    async_foo(2)
    async_foo(3)
end

async_foo(1)

entry()
bar(1 ,2)

function test1()
    @sync for i in 1:10
        Threads.@spawn bar(i, 1)
    end
end
function test2()
    Threads.@threads for i in 1:10
        bar(i, 1)
    end
end
@time test1()
@time test2()

# just sync
@sync begin end
quote
    let var"##sync#41" = Base.Channel(Base.Inf)
        var"#862#v" = begin
        end
        Base.sync_end(var"##sync#41")
        var"#862#v"
    end
end

# just spawn
Threads.@spawn bar(i, 1)
quote
    let
        local _task = Base.Threads.Task((()->begin
                                                      bar(i, 1)
                                                  end))
        (_task).sticky = false
        if $(Expr(:islocal, Symbol("##sync#41")))
            Base.Threads.put!(var"##sync#41", _task)
        end
        Base.Threads.schedule(_task)
        _task
    end
end

# sync spawn
quote
    begin
        let var"##sync#41" = Base.Channel(Base.Inf)
            _v = for i = 1:10
                begin
                    let
                        local _task = Base.Threads.Task((()->begin
                                                                      bar(i, 1)
                                                                  end))
                        (_task).sticky = false
                        if $(Expr(:islocal, Symbol("##sync#41")))
                            Base.Threads.put!(var"##sync#41", _task)
                        end
                        Base.Threads.schedule(_task)
                        _task
                    end
                end
            end
            Base.sync_end(var"##sync#41")
            _v
        end
    end
end

# thread for
quote
    begin
        local _threadsfor_fun
        let _range = 1:10
            function _threadsfor_fun(_onethread = false)
                _r = _range
                _lenr = Base.Threads.length(_r)
                if _onethread
                    _tid = 1
                    (_len, _rem) = (_lenr, 0)
                else
                    _tid = Base.Threads.threadid()
                    (_len, _rem) = Base.Threads.divrem(_lenr, Base.Threads.nthreads())
                end
                if _len == 0
                    if _tid > _rem
                        return
                    end
                    (_len, _rem) = (1, 0)
                end
                _f = Base.Threads.firstindex(_r) + (_tid - 1) * _len
                _l = (_f + _len) - 1
                if _rem > 0
                    if _tid <= _rem
                        _f = _f + (_tid - 1)
                        _l = _l + _tid
                    else
                        _f = _f + _rem
                        _l = _l + _rem
                    end
                end
                for _i = _f:_l
                    local i = begin
                        $(Expr(:inbounds, true))
                        local _val = _r[_i]
                        $(Expr(:inbounds, :pop))
                        _val
                    end
                    begin
                        bar(i, 1)
                    end
                end
            end
        end
        if Base.Threads.threadid() != 1 || ccall(:jl_in_threaded_region, Base.Threads.Cint, ()) != 0
            (Base.Threads.Base).invokelatest(_threadsfor_fun, true)
        else
            Base.Threads.threading_run(_threadsfor_fun)
        end
        Base.Threads.nothing
    end
end
