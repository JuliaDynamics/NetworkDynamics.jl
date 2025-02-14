module NetworkDynamicsInspector

using Bonito
using Bonito.Observables
using Bonito: Grid, @js_str, onload, jsrender

export ContinuousSlider
include("widgets.jl")

# function inspect(sol)
#     t = Observable(0.0)

# end

# function gparguments(sol;
#                      t,
#                      state,
#                      ncolorrange,
#                      ncolorscheme,
#                      n_rel_to_u0,
#                      sel_nodes = Observable(Set{Int}()))

#     t::Observable = t isa Observable ? t : Observable(t)
#     NV = nv(network)

#     node_marker = try
#         markdict = Dict(:load => :rect, :gen => :circle, :syncon => :circle);
#         [markdict[k] for k in get_prop(network, 1:NV, :type)];
#     catch
#         markerset = [:circle, :rect, :utriangle, :cross, :diamond, :dtriangle, :pentagon, :xcross]
#         groups = sol.prob.f.f.unique_v_indices
#         markers = Vector{Symbol}(undef, nv(network))
#         for (i,g) in enumerate(groups)
#             markers[g] .= markerset[i]
#         end
#         markers
#     end

#     u0statevec = Observable(Vector{Float64}(undef, NV))

#     on(nstatelens; update=true) do lens
#         u0statevec[] .= @views lens(sol.t)[1, :]
#     end

#     statevec = Observable(Vector{Float64}(undef, NV))

#     onany(t, nstatelens, n_rel_to_u0) do t, lens, rel
#         statevec[] .= @views lens(t)[1, :]
#         if rel
#             statevec[] .-= u0statevec[]
#         end
#         notify(statevec)
#     end

#     node_color = Observable(Vector{RGB{Float64}}(undef, NV))
#     onany(statevec, ncolorrange, ncolorscheme) do statevec, range, scheme
#         for i in 1:NV
#             node_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
#         end
#         notify(node_color)
#     end
#     notify(t)

#     SMALL = 30
#     BIG = 50
#     node_size = Observable(fill(SMALL, NV))
#     onany(sel_nodes) do selected
#         fill!(node_size[], SMALL)
#         for sel in selected
#             node_size[][sel] = BIG
#         end
#         notify(node_size)
#     end
#     notify(sel_nodes)

#     return (;layout=read_pos_or_spring,
#             node_marker,
#             node_color,
#             node_size,
#             node_attr=(;colorrange=ncolorrange, colormap=ncolorscheme));
# end

end # module NetworkDynamicsInspector
