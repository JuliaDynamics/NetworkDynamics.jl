using OrderedCollections
using AbstractTrees
using Chairmarks
using PrettyTables

struct BenchmarkDict{D}
    d::D
end

struct SampleComparison
    target::Chairmarks.Sample
    baseline::Chairmarks.Sample
end

hastarget(::SampleComparison) = true
hasbaseline(::SampleComparison) = true
hastarget(::Chairmarks.Sample) = true
hasbaseline(::Chairmarks.Sample) = false
gettarget(sc::SampleComparison) = sc.target
getbaseline(sc::SampleComparison) = sc.baseline
gettarget(s::Chairmarks.Sample) = s

BenchmarkDict() = BenchmarkDict(OrderedDict())
BenchmarkDict(args...) = BenchmarkDict(OrderedDict(args...))

Base.isempty(bd::BenchmarkDict) = isempty(bd.d)
Base.keys(bd::BenchmarkDict) = keys(bd.d)
Base.values(bd::BenchmarkDict, k) = values(bd.d)
Base.getindex(bd::BenchmarkDict, k) = get!(bd.d, k, BenchmarkDict())
Base.getindex(bd::BenchmarkDict, k...) = bd[k[1:(end - 1)]...][k[end]]
Base.setindex!(bd::BenchmarkDict, v, k) = setindex!(bd.d, v, k)
Base.setindex!(bd::BenchmarkDict, v, k...) = setindex!(bd[k[1:(end - 1)]...], v, k[end])

AbstractTrees.children(bd::BenchmarkDict) = pairs(bd.d)
AbstractTrees.printnode(io::IO, bd::BenchmarkDict) = print(io, "BenchmarkDict()")

Base.show(io::IO, ::MIME"text/plain", bd::BenchmarkDict) = pretty_table(io, bd)

function PrettyTables.pretty_table(io::IO, bd::BenchmarkDict; kwargs...)
    flatk, flatv = _flatten(bd)
    if isempty(flatk)
        print(io, "BenchmarkDict()")
        return
    end

    keycols = [[x[i] for x in flatk] for i in eachindex(flatk[1])]
    for col in keycols
        for i in length(col):-1:2
            if col[i] == col[i-1]
                col[i] = ""
            end
        end
    end
    ttime = map(flatv) do _v
        v = gettarget(_v)
        @assert v isa Chairmarks.Sample
        buf = IOBuffer()
        Chairmarks.print_time(buf, v.time)
        String(take!(buf))
    end
    tallocs = map(flatv) do _v
        v = gettarget(_v)
        @assert v isa Chairmarks.Sample
        v.allocs > 0 ? repr(v.allocs) : ""
    end

    data = hcat(keycols..., ttime, tallocs)
    header = vcat("Key", ["" for i in 1:length(keycols)-1]..., "Time", "Allocs")

    if any(hasbaseline, flatv)
        btime = map(flatv) do _v
            hasbaseline(_v) || return ""
            v = getbaseline(_v)
            @assert v isa Chairmarks.Sample
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
            @assert v isa Chairmarks.Sample
            v.allocs > 0 ? repr(v.allocs) : ""
        end
        callocs = map(flatv) do v
            hasbaseline(v) || return ""
            gettarget(v).allocs == getbaseline(v).allocs && return 0
            Δ = (gettarget(v).allocs - getbaseline(v).allocs) / getbaseline(v).allocs * 100
            # round(Δ,digits=1)
        end
        data = hcat(keycols..., ttime, btime, ctime, tallocs, ballocs, callocs)
        header = (vcat("Key", ["" for i in 1:length(keycols)-1]..., "Time", "", "","Allocs","",""),
                  vcat(["" for i in 1:length(keycols)]..., "target", "baseline", "","target","baseline",""))
        hl_bad = Highlighter((data, i, j) -> (j ∈ length(keycols) .+ [3,6]) && (data[i, j] > 0), crayon"red bold");
        hl_good = Highlighter((data, i, j) -> (j ∈ length(keycols) .+ [3,6]) && (data[i, j] < 0), crayon"green bold");
        kwargs = (kwargs..., highlighters=(hl_bad, hl_good), formatters    = ft_printf("%+5.1f %%", length(keycols) .+ [3,6]))
    end

    pretty_table(io, data; header, header_alignment=:l, kwargs...)
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

function compare(target, baseline)
    tkeys, _ = _flatten(target)
    bkeys, _ = _flatten(baseline)

    inters = tkeys ∩ bkeys
    bd = BenchmarkDict()
    for k in inters
        bd[k...] = SampleComparison(target[k...], baseline[k...])
    end
    bd
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

target = deserialize("./2024-06-05_124146_target.data")
baseline = deserialize("./2024-06-05_124146_baseline.data")
compare(target, baseline)


end
