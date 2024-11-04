using OrderedCollections
using AbstractTrees
using Chairmarks
using PrettyTables
using NetworkDynamics
using Test
using Printf

struct BenchmarkDict{D}
    d::D
end

struct BenchmarkResult
    time::Float64
    allocs::Int
    samples::Int
    value::Any
end
function BenchmarkResult(be::Chairmarks.Benchmark, value=nothing)
    res = minimum(be)
    BenchmarkResult(res.time, res.allocs, length(be.samples), value)
end

struct SampleComparison
    target::BenchmarkResult
    baseline::BenchmarkResult
end

hastarget(::SampleComparison) = true
hasbaseline(::SampleComparison) = true
hastarget(::BenchmarkResult) = true
hasbaseline(::BenchmarkResult) = false
gettarget(sc::SampleComparison) = sc.target
getbaseline(sc::SampleComparison) = sc.baseline
gettarget(s::BenchmarkResult) = s

BenchmarkDict() = BenchmarkDict(OrderedDict())
BenchmarkDict(args...) = BenchmarkDict(OrderedDict(args...))

Base.isempty(bd::BenchmarkDict) = isempty(bd.d)
Base.keys(bd::BenchmarkDict) = Base.keys(bd.d)
Base.values(bd::BenchmarkDict) = Base.values(bd.d)
Base.getindex(bd::BenchmarkDict, k) = get!(bd.d, k, BenchmarkDict())
Base.getindex(bd::BenchmarkDict, k...) = bd[k[1:(end - 1)]...][k[end]]
Base.setindex!(bd::BenchmarkDict, v, k) = setindex!(bd.d, v, k)
Base.setindex!(bd::BenchmarkDict, v, k...) = setindex!(bd[k[1:(end - 1)]...], v, k[end])
Base.haskey(bd::BenchmarkDict, k) = haskey(bd.d, k)
Base.haskey(bd::BenchmarkDict, k...) = haskey(bd.d, k[1]) && haskey(bd[k[1]], k[2:end]...)

AbstractTrees.children(bd::BenchmarkDict) = pairs(bd.d)
AbstractTrees.printnode(io::IO, bd::BenchmarkDict) = print(io, "BenchmarkDict()")

Base.show(io::IO, ::MIME"text/plain", bd::BenchmarkDict) = pretty_table(io, bd)

function PrettyTables.pretty_table(io::IO, bd::BenchmarkDict; kwargs...)
    flatk, flatv = _flatten(bd)
    if isempty(flatk)
        print(io, "BenchmarkDict()")
        return
    end

    numbers = map(flatk) do line
        line[end]
    end

    keycols = if numbers isa AbstractVector{<:Number}
        basekey = map(flatk) do line
            ret = ""
            for kp in line[1:end-2]
                ret *= kp*" → "
            end
            ret *= line[end-1]
        end
        basekey, numbers
    else
        [[i∈eachindex(x) ? x[i] : "" for x in flatk] for i in 1:maximum(length.(flatk))]
    end

    for col in keycols
        for i in length(col):-1:2
            if col[i] == col[i-1]
                col[i] = ""
            end
        end
    end
    ttime = map(flatv) do _v
        v = gettarget(_v)
        @assert v isa BenchmarkResult
        buf = IOBuffer()
        Chairmarks.print_time(buf, v.time)
        String(take!(buf))
    end
    tallocs = map(flatv) do _v
        v = gettarget(_v)
        @assert v isa BenchmarkResult
        v.allocs > 0 ? repr(v.allocs) : ""
    end
    samples = map(flatv) do _v
        v = gettarget(_v)
        @assert v isa BenchmarkResult
        v.samples
    end

    data = hcat(keycols..., ttime, tallocs, samples)
    header = vcat("Key", ["" for i in 1:length(keycols)-1]..., "Time", "Allocs", "Samples")

    if any(hasbaseline, flatv)
        btime = map(flatv) do _v
            hasbaseline(_v) || return ""
            v = getbaseline(_v)
            @assert v isa BenchmarkResult
            buf = IOBuffer()
            Chairmarks.print_time(buf, v.time)
            String(take!(buf))
        end
        ctime = map(flatv) do v
            hasbaseline(v) || return ""
            Δ = (gettarget(v).time - getbaseline(v).time) / getbaseline(v).time * 100
            # round(Δ,digits=1)
        end
        ballocs = map(flatv) do _v
            hasbaseline(_v) || return ""
            v = getbaseline(_v)
            @assert v isa BenchmarkResult
            v.allocs > 0 ? repr(v.allocs) : ""
        end
        callocs = map(flatv) do v
            hasbaseline(v) || return ""
            gettarget(v).allocs == getbaseline(v).allocs && return 0
            Δ = (gettarget(v).allocs - getbaseline(v).allocs) / getbaseline(v).allocs * 100
            # round(Δ,digits=1)
        end

        _kwdict = Dict(kwargs)
        if haskey(_kwdict, :backend) && _kwdict[:backend] == Val(:markdown)
            ctime = map(ctime) do num
                sym = if num ≤ 0
                    "✅"
                elseif num > 5
                    "❌"
                else
                    "✔"
                end
                repr(round(num, digits=2)) * " % " * sym
            end
            callocs = map(callocs) do num
                sym = if num ≤ 0
                    "✅"
                elseif num > 5
                    "❌"
                else
                    "✔"
                end
                repr(round(num, digits=2)) * " % " * sym
            end
        else
            hl_bad = Highlighter(crayon"red bold") do data, i, j
                (j ∈ length(keycols) .+ [3,6]) && data[i, j] isa Number && data[i,j] > 0
            end
            hl_good = Highlighter(crayon"green bold") do data, i, j
                (j ∈ length(keycols) .+ [3,6]) && data[i, j] isa Number && (data[i, j] < 0)
            end
            formatters = (v, _, col) -> begin
                if col in length(keycols) .+ [3,6]
                    if v < -5
                        @sprintf "%+5.1f %% ✅" v
                    elseif v < 5
                        @sprintf "%+5.1f %% ➖" v
                    else
                        @sprintf "%+5.1f %% ❌" v
                    end
                else
                    v
                end
            end
            kwargs = (; kwargs..., highlighters=(hl_bad, hl_good), formatters)
        end

        data = hcat(keycols..., ttime, btime, ctime, tallocs, ballocs, callocs)
        header = (vcat("Key", ["" for i in 1:length(keycols)-1]..., "Time", "", "","Allocs","",""),
                  vcat(["" for i in 1:length(keycols)]..., "target", "baseline", "","target","baseline",""))
    end

    alignment = vcat([:l for _ in 1:length(keycols)-1]..., [:r for _ in 1:size(data,2)-length(keycols)+1]...)
    pretty_table(io, data; header, header_alignment=:l, alignment, kwargs...)

end

function test_return_values(bd)
    @testset "Compare Results" begin
        flatk, flatv = _flatten(bd)
        for (k,b) in zip(flatk, flatv)
            b isa BenchmarkResult && continue
            bsv = getbaseline(b).value
            tgv = gettarget(b).value
            @test bsv == tgv || isapprox(bsv, tgv)
            if !(bsv == tgv || isapprox(bsv, tgv))
                @warn "Values differ for $(k)! $(extrema(bsv-tgv))"
            end
        end
    end
end

function _flatten(bd::BenchmarkDict, prekey=[], _flatk=[], _flatv=[])
    for (k, v) in bd.d
        if v isa BenchmarkDict
            _flatten(v, [prekey; k], _flatk, _flatv)
        else
            _flatk = push!(_flatk, [prekey; k])
            _flatv = push!(_flatv, v)
        end
    end
    _flatk, _flatv
end

function compare(target, baseline; alltarget=false)
    tkeys, _ = _flatten(target)
    bd = BenchmarkDict()

    for k in tkeys
        if haskey(baseline, k...)
            bd[k...] = SampleComparison(target[k...], baseline[k...])
        elseif alltarget
            bd[k...] = target[k...]
        end
    end
    bd
end

function plot_over_N(target, baseline=nothing)
    # exs = ["seq_buf", "ka_buf, poly_buf, threaded_buf", "seq", "ka", "poly", "threaded"]
    # exs = ["seq_buf", "poly_buf", "threaded_buf", "ka_buf"]
    # exs = ["seq", "seq_buf"]
    # exs = ["ka", "ka_buf"]
    # exs = ["poly", "poly_buf"]
    # exs = ["threaded", "threaded_buf"]
    # exs = ["seq"]
    # aggrs = ["nnlib","KA","seq","poly","thrd"]
    # aggrs = ["seq"]

    if baseline != nothing
        comp = compare(target, baseline; alltarget=true)
    else
        comp = target
    end

    fig = Makie.Figure(size=(2000,2000))

    bmkeys = [
        ("diffusion", "static_edge"),
        ("diffusion", "ode_edge"),
        ("kuramoto", "homogeneous"),
        ("kuramoto", "heterogeneous"),
    ]

    for (row, key) in pairs(bmkeys)
        dat = comp[key...]

        allex = sort(collect(filter(!isequal("assemble"), keys(dat))))
        defex = "poly_buf" ∈ allex ? "poly_buf" : allex[1]
        allagg = sort(collect(keys(dat[defex])))
        defagg = "seq" ∈ allagg ? "seq" : allagg[1]

        ax = Makie.Axis(fig[row,1]; xscale=log10, yscale=log10, ylabel="coreloop time", title="$key executions agg=$defagg")
        for ex in allex
            datagg = dat[ex, defagg]
            isempty(datagg) && continue
            N = collect(keys(datagg))
            ttime = getproperty.(gettarget.(values(datagg)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex")
            if all(hasbaseline.(values(datagg)))
                btime = getproperty.(getbaseline.(values(datagg)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
        # gpu timeseries
        for (ex, agg) in [("ka_buf_cuda32", "sprs_cuda32"), ("ka_buf_cuda64", "sprs_cuda64")]
            datagg = dat[ex, agg]
            isempty(datagg) && continue
            N = collect(keys(datagg))
            ttime = getproperty.(gettarget.(values(datagg)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex")
            if all(hasbaseline.(values(datagg)))
                btime = getproperty.(getbaseline.(values(datagg)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end

        Makie.axislegend(ax; position=:lt)

        ax = Makie.Axis(fig[row,2]; xscale=log10, yscale=log10, ylabel="coreloop time", title="$key aggregations ex=$defex")
        for agg in allagg
            datagg = dat[defex, agg]
            isempty(datagg) && continue
            N = collect(keys(datagg))
            ttime = getproperty.(gettarget.(values(datagg)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$agg")
            if all(hasbaseline.(values(datagg)))
                btime = getproperty.(getbaseline.(values(datagg)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
        Makie.axislegend(ax; position=:lt)
    end

    fig
end

# code for interactive testing
@static if false
bd = BenchmarkDict()
bd["Group 1", "Subgroup 1", "Benchmark 1"] = BenchmarkResult(1.0, 100)
bd["Group 1", "Subgroup 1", "Benchmark 2"] = BenchmarkResult(2.0, 100)
bd["Group 1", "Subgroup 2", "Benchmark 1"] = BenchmarkResult(3.0, 100)
bd["Group 1", "Subgroup 2", "Benchmark 2"] = BenchmarkResult(4.0, 100)
bd["Group 2", "Subgroup 1", "Benchmark 1"] = BenchmarkResult(5.0, 100)
bd["Group 2", "Subgroup 1", "Benchmark 2"] = BenchmarkResult(6.0, 100)
bd["Group 2", "Subgroup 2", "Benchmark 1"] = BenchmarkResult(7.0, 100)
bd["Group 2", "Subgroup 2", "Benchmark 2"] = BenchmarkResult(9.0, 100)

target = deserialize(sort(filter(contains("target.data"), readdir()))[end])
baseline = deserialize(sort(filter(contains("baseline.data"), readdir()))[end])
comp = compare(target, baseline)
test_return_values(comp)

entry =  comp["diffusion","ode_edge","ka_buf","KA",6]
entry =  comp["diffusion","ode_edge","seq_buf","seq",6]
entry.target.value[1] - entry.baseline.value[1]
entry.target.value[2]
entry.baseline.value[2]

entry =  comp["kuramoto","heterogeneous","ka_buf","KA",6]
# entry =  comp["kuramoto","heterogeneous","seq_buf","seq",6]
entry.target.value[1]
entry.baseline.value[1]
entry.target.value[2]
entry.baseline.value[2]

target = deserialize(sort(filter(contains("target.data"), readdir()))[end])
baseline = deserialize(sort(filter(contains("baseline.data"), readdir()))[end])
comp = compare(target,baseline)

pretty_table(comp["diffusion","static_edge"]; backend=Val(:markdown))

end
