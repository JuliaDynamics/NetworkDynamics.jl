# Chapter 1 - Pure Julia
begin
    using Graphs, OrdinaryDiffEq

    N = 10
    g = watts_strogatz(N, 2, 0.0)
    const B = incidence_matrix(g; oriented=true)
    const B_t = transpose(B)

    function kuramoto_network!(dθ, θ, ω, t)
        dθ .= ω .- 5.0 .* (B * sin.(B_t * θ))
        return nothing
    end

    ω = (collect(1:N) .- sum(1:N) / N) / N
    x0 = (collect(1:N) .- sum(1:N) / N) / N
    tspan = (0.0, 4.0)
    prob = ODEProblem(kuramoto_network!, x0, tspan, ω)

    sol = solve(prob, Tsit5())
end

# Chapter 2 - Homogeneous ND

begin
    function kuramoto_edge!(edge, θ_s, θ_d, σ, t)
        edge .= σ .* sin.(θ_s .- θ_d)
        return nothing
    end
    function kuramoto_vertex!(dθ, θ, edges, ω, t)
        dθ .= ω
        for edge in edges
            dθ[1] += edge[1]
        end
        return nothing
    end


    using NetworkDynamics

    f_node! = ODEVertex(; f=kuramoto_vertex!, dim=1, sym=[:θ])
    f_edge! = StaticEdge(; f=kuramoto_edge!, dim=1)
    nd! = network_dynamics(f_node!, f_edge!, g)


    vertexp = ω
    edgep   = 5.0
    p       = (vertexp, edgep)

    nd_prob = ODEProblem(nd!, x0, tspan, p)
    nd_sol = solve(nd_prob, Tsit5())


    using Plots

    plot(nd_sol; ylabel="θ")
end
savefig("kuramoto_homogeneous.pdf")



# Hidden Color configuration
begin
    membership = Int64.(ones(N))
    membership[1] = 2
    membership[N÷2] = 3
    nodecolor = [colorant"lightseagreen", colorant"orange", colorant"darkred"]
    nodefillc = nodecolor[membership]
    nodefillc = reshape(nodefillc, 1, N)
end

# Chapter 3 - Heterogeneous ND
begin
    function kuramoto_inertia!(dv, v, edges, p, t)
        dv[1] = v[2]
        dv[2] = p - v[2]
        for edge in edges
            dv[2] += edge[1]
        end
        return nothing
    end

    inertia! = ODEVertex(; f=kuramoto_inertia!, dim=2, sym=[:θ, :ω])

    static! = StaticVertex(; f=(θ, edges, c, t) -> θ .= c,
                           dim=1, sym=[:θ])

    function kuramoto_edge!(edge, v_s, v_d, σ, t)
        edge[1] = σ * sin(v_s[1] - v_d[1])
        return nothing
    end

    v_arr = Array{VertexFunction}([f_node! for v in vertices(g)])
    v_arr[1] = inertia!
    v_arr[N÷2] = static!
    nd_hetero! = network_dynamics(v_arr, f_edge!, g)

    # Parameters and inital conditions

    insert!(x0, 2, 3.0)
    prob_hetero = ODEProblem(nd_hetero!, x0, tspan, p)
    sol_hetero = solve(prob_hetero, Rodas4())

    vars = syms_containing(nd_hetero!, :θ)
    plot(sol_hetero; ylabel="θ", vars=vars, lc=nodefillc)
end
savefig("kuramoto_heterogeneous.pdf")


# Chapter 4 - Fancy ND (Delays)


begin
    function kuramoto_delay_edge!(edge, v_s, v_d, h_v_s, h_v_d, p, t)
        edge[1] = p * sin(v_s[1] - h_v_d[1])
        return nothing
    end
    dedge! = StaticDelayEdge(; f=kuramoto_delay_edge!, dim=1)

    using DelayDiffEq

    h(out, p, t) = (out .= x0)
    τ = 0.1
    p = (vertexp, edgep, τ)
    nd_delay! = network_dynamics(v_arr, dedge!, g)
    prob_delay = DDEProblem(nd_delay!, x0, h, tspan, p; constant_lags=[τ])

    sol_delay = solve(prob_delay, MethodOfSteps(Rodas4(; autodiff=false)))

    plot(sol_delay; ylabel="θ", vars=vars, lc=nodefillc)
end

savefig("kuramoto_delay.pdf")


# Chapter  5 - Fancy ND II (Callbacks)

begin
    using DiffEqCallbacks

    θ_idxs = idx_containing(nd_delay!, :θ)

    function condition(out, u, t, integrator)
        out .= (u[θ_idxs] .- 0.5) .*
               (u[θ_idxs] .+ 0.5)
        return nothing
    end

    function affect!(integrator, idx)
        stable_edges = map(edg -> idx ∉ edg, Pair.(edges(g)))
        integrator.p = (integrator.p[1], stable_edges .* integrator.p[2], integrator.p[3])
        return nothing
    end

    cb = VectorContinuousCallback(condition, affect!, 10)
    prob_cb = remake(prob_delay;
                     p=(vertexp, edgep .* ones(N), τ))
    sol_cb = solve(prob_cb, MethodOfSteps(Rodas4(; autodiff=false)); callback=cb)



    plot(sol_cb; ylabel="θ", vars=vars, lc=nodefillc)
    hline!([-.5]; color=:black, width=1.0, line=:dot, label="")
    hline!([0.5]; color=:black, width=1.0, line=:dot, label="")
end

savefig("kuramoto_callback.pdf")
