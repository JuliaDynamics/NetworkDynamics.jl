"""
For a given model, collect all referenced models
"""
collect_model_dependencies(model) = unique!(_collect_model_dependencies(String[], model))
_collect_model_dependencies(deps, other) = deps
function _collect_model_dependencies(deps, model::Union{OrderedDict,Vector})
    for val in values(model)
        if val isa String
            m = match(r"^Models.(.*)$", val)
            if !isnothing(m)
                push!(deps, m[1])
            end
        else
            _collect_model_dependencies(deps, val)
        end
    end
    deps
end

"""
    recursiv_resolve!(models, container)

Go through container, and replace every occurence of "Models.<modelname>" with
deepcopy of actual model dict.
"""
recursive_resolve!(models, _) = nothing
function recursive_resolve!(models, container::Union{OrderedDict,Vector})
    for (key, value) in pairs(container)
        if value isa String
            m = match(r"^Models.(.*)$", value)
            if !isnothing(m)
                container[key] = deepcopy(models[m[1]])
            end
        else
            recursive_resolve!(models, value)
        end
    end
    nothing
end

"""
    update_parent_property!(parent, key, value)

Update the property in the parent structure.
- if key is nested (i.e. "a.b.c"), it will update parent["a"]["b"]["c"]
- if key is indexed (i.e. "a[1]"), it will update parent["a"][1]
"""
function update_parent_property!(parent, key, value)
    @info "Update $key -> $value"
    if contains(key, ".")
        key_head, key_tail = split(key, ".", limit=2)
        @info "Found nested key $key_head -> $key_tail"
        m = match(r"^(.*)\[(.*)\]$", key_head)
        _parent = if !isnothing(m)
            parent[m[1]][parse(Int, m[2])]
        else
            parent[key]
        end
        @info m
        update_parent_property!(_parent, key_tail, value)
    else
        m = match(r"^(.*)\[(.*)\]$", key)
        if !isnothing(m)
            if !(value isa eltype(parent[m[1]]))
                # need to widen type
                parent[m[1]] = collect(Any, parent[m[1]])
            end
            parent[m[1]][parse(Int, m[2])] = value
        else
            parent[key] = value
        end
    end
end
"""
    merge_properties!(parent, child)

Takes a child dict, and applies every key-value pair to the parent dict.
keys will be resolved using update_parent_property!
"""
function merge_properties!(parent, child)
    for (key, value) in pairs(child)
        update_parent_property!(parent, key, value)
    end
    parent
end

"""
    depth_first_flatten!(container)

This goes through the container, and if it finds a dict with key "MODEL", it will
merge the properties of the referenced parent model with the current dict.
"""
depth_first_flatten!(_) = nothing
function depth_first_flatten!(container::Union{OrderedDict,Vector})
    for (key, sub) in pairs(container)
        depth_first_flatten!(sub)
        if sub isa Union{OrderedDict} && haskey(sub, "MODEL")
            parent_model = sub["MODEL"]
            delete!(sub, "MODEL")
            container[key] = merge_properties!(parent_model, sub)
        end
    end
    nothing
end

function parse_network(file)
    data = YAML.load_file(file, dicttype=OrderedDict{Any,Any})
    dependencies = OrderedDict(keys(data["Models"]) .=> map(collect_model_dependencies, values(data["Models"])))

    # dependencies muss by subset of models
    @assert reduce(vcat, values(dependencies)) ⊆ keys(data["Models"])
    resolved_models = findall(deps -> isempty(deps), dependencies)
    delete!.(Ref(dependencies), resolved_models) # delete resolve dependnecies

    # resolve all models, i.e. replace parent/template models with full dict
    while true
        isempty(dependencies) && break
        resolved_something = false
        for (mod, deps) in dependencies
            if deps ⊆ resolved_models
                @info "Resolving $mod -> $deps"
                recursive_resolve!(data["Models"], data["Models"][mod])
                delete!(dependencies, mod)
                push!(resolved_models, mod)
                resolved_something = true
            else
                @info "Skip $mod -> $deps"
            end
        end
        resolved_something && continue
        error("Could not resolve all dependencies")
    end

    # with all deps resolve, we can replace the models edge and vertex list
    recursive_resolve!(data["Models"], data["Vertices"])
    recursive_resolve!(data["Models"], data["Edges"])

    # traverses the tree and flattens the model properties, i.e. merges
    # the local given properties with the parent/template model properties
    depth_first_flatten!(data)
    return data
end

recursive_construct(constructors, data) = _depth_first_construct!(constructors, deepcopy(data))
_depth_first_construct!(_, data) = data
_depth_first_construct!(_, s::String) = startswith(s, ":") ? Symbol(s[2:end]) : s
function _depth_first_construct!(constructors, data::Vector)
    if data isa Vector && eltype(data) !== Any
        data = convert(Vector{Any}, data)
    end
    for (key, value) in pairs(data)
        data[key] = _depth_first_construct!(constructors, value)
    end
    data
end
function _depth_first_construct!(constructors, data::OrderedDict)
    if data isa OrderedDict && eltype(values(data)) !== Any
        data = convert(OrderedDict{eltype(keys(data)), Any}, data)
    end

    for (key, value) in pairs(data)
        data[key] = _depth_first_construct!(constructors, value)
    end
    if haskey(data, "CONSTRUCTOR")
        if !haskey(constructors, data["CONSTRUCTOR"])
            error("No constructor found for $(data["CONSTRUCTOR"])")
        end
        args = get(data, "ARGS", [])
        kwargs = Dict{Symbol, Any}()
        for (key, value) in data
            key ∈ ("CONSTRUCTOR", "ARGS") && continue
            kwargs[Symbol(key)] = value
        end
        return constructors[data["CONSTRUCTOR"]](args...; kwargs...)
    end
    data
end

function build_network(data, constructors)
    vertexm = VertexModel[]
    for (k, v) in data["Vertices"]
        vm = recursive_construct(constructors, v)
        set_graphelement!(vm, k)
        push!(vertexm, vm)
    end
    edgem   = EdgeModel[]
    for (k, v) in data["Vertices"]
        vm = recursive_construct(constructors, v)

        m = match(r"^(.*)=>(.*)$", k)
        isnothing(m) && error("Edge key must be of form 'src=>dst', got $k")
        src = tryparse(Int, m[1])
        dst = tryparse(Int, m[2])
        isnothing(src) && src = m[1]
        isnothing(dst) && src = m[1]
        set_graphelement!(vm, src, dst)
        push!(vertexm, vm)
    end
end
