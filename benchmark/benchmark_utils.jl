using OrderedCollections
using AbstractTrees
using Chairmarks
using PrettyTables
using Test

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
            formatters = ft_printf("%+5.1f %%", length(keycols) .+ [3,6])
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

function plot_over_N(target, baseline)
    exs = ["seq_buf", "ka_buf"]
    aggrs = #=["nnlib",=#["KA","seq","poly"]

    comp = compare(target, baseline; alltarget=true)

    fig = Makie.Figure(size=(2000,1000))
    ax = Makie.Axis(fig[1,1]; xscale=log10, yscale=log10, ylabel="coreloop time", title="Diffusion - static_edge")

    for ex in exs
        for aggr in aggrs
            dat = comp["diffusion", "static_edge", ex, aggr]
            N = collect(keys(dat))
            ttime = getproperty.(gettarget.(values(dat)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex $aggr")
            if all(hasbaseline.(values(dat)))
                btime = getproperty.(getbaseline.(values(dat)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
    end
    Makie.axislegend(ax; position=:lt)

    ax = Makie.Axis(fig[1,2]; xscale=log10, yscale=log10, ylabel="coreloop time", title="Diffusion - ode_edge")
    for ex in exs
        for aggr in aggrs
            dat = comp["diffusion", "ode_edge", ex, aggr]
            N = collect(keys(dat))
            ttime = getproperty.(gettarget.(values(dat)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex $aggr ")
            if all(hasbaseline.(values(dat)))
                btime = getproperty.(getbaseline.(values(dat)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
    end
    Makie.axislegend(ax; position=:lt)

    ax = Makie.Axis(fig[2,1]; xscale=log10, yscale=log10, ylabel="coreloop time", title="kuramoto - homogeneous")
    for ex in exs
        for aggr in aggrs
            dat = comp["kuramoto", "homogeneous", ex, aggr]
            N = collect(keys(dat))
            ttime = getproperty.(gettarget.(values(dat)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex $aggr ")
            if all(hasbaseline.(values(dat)))
                btime = getproperty.(getbaseline.(values(dat)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
    end
    Makie.axislegend(ax; position=:lt)

    ax = Makie.Axis(fig[2,2]; xscale=log10, yscale=log10, ylabel="coreloop time", title="kuramoto - heterogeneous")
    for ex in exs
        for aggr in aggrs
            dat = comp["kuramoto", "heterogeneous", ex, aggr]
            N = collect(keys(dat))
            ttime = getproperty.(gettarget.(values(dat)), :time)
            sc = Makie.scatterlines!(ax, N, ttime, label="$ex $aggr ")
            if all(hasbaseline.(values(dat)))
                btime = getproperty.(getbaseline.(values(dat)), :time)
                Makie.scatterlines!(ax, N, btime; linestyle=:dash, color=sc.color)
            end
        end
    end
    Makie.axislegend(ax; position=:lt)
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
