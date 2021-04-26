module nd_ODE_Static_mod

using ..NetworkStructures
using ..ComponentFunctions
using ..Utilities


export nd_ODE_Static
export StaticEdgeFunction
export prep_gd


# If the type of dx and x match, we swap out v_array for x
@inline function prep_gd(dx::AbstractArray{T}, x::AbstractArray{T}, gd::GraphData{GDB, T, T}, gs) where {GDB, T}
    # We construct views into an Array of size dim_v, so when we swap the
    # underlying Array it has to have to same size
    if size(x) == (gs.dim_v,)
         swap_v_array!(gd, x)
         return gd
    else
         error("Size of x does not match the dimension of the system.")
    end
end

# If the type of dx and x do not match, we swap initialize a new GraphData object
# that is based on the type of dx for the edge buffer.
# Some solvers take the derivative with respect to time, thus x will not be dual
# but dx will be, leading to errors otherwise
@inline function prep_gd(dx, x, gd, gs)
    # Type mismatch
    if size(x) == (gs.dim_v,)
        e_array = similar(dx, gs.dim_e)
        return GraphData(x, e_array, gs)
    else
        error("Size of x does not match the dimension of the system.")
   end
end



@Base.kwdef struct nd_ODE_Static{G, GD, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G #redundant?
    graph_structure::GraphStruct
    graph_data::GD
    parallel::Bool # enables multithreading for the core loop
end


function (d::nd_ODE_Static)(dx, x, p, t)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)
    gd = prep_gd(dx, x, d.graph_data, d.graph_structure)

    @nd_threads d.parallel for i in 1:gs.num_e
        maybe_idx(d.edges!, i).f!(
            get_edge(gd, i),
            get_src_vertex(gd, i),
            get_dst_vertex(gd, i),
            p_e_idx(p, i),
            t)
    end

    @assert size(dx) == size(x) "Sizes of dx and x do not match"

    @nd_threads d.parallel for i in 1:gs.num_v
        maybe_idx(d.vertices!,i).f!(
            view(dx,gs.v_idx[i]),
            get_vertex(gd, i),
            get_dst_edges(gd, i),
            p_v_idx(p, i),
            t)
    end

    nothing
end


function (d::nd_ODE_Static)(x, p, t, ::Type{GetGD})
    gd = prep_gd(x, x, d.graph_data, d.graph_structure)
    gs = d.graph_structure
    checkbounds_p(p, gs.num_v, gs.num_e)

    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
        maybe_idx(d.edges!,i).f!(
            get_edge(gd, i),
            get_src_vertex(gd, i),
            get_dst_vertex(gd, i),
            p_e_idx(p, i),
            t)
    end

    gd
end


function (d::nd_ODE_Static)(::Type{GetGS})
    d.graph_structure
end


# For compatibility with PowerDynamics

struct StaticEdgeFunction
    nd_ODE_Static
end

function (sef::StaticEdgeFunction)(x, p, t)
    gd = sef.nd_ODE_Static(x, p, t, GetGD)

    # Allocating but non-breaking. No special knowlage about the fields of
    # GraphData needed anymore. Is there a real need for this in PD.jl?
    num_v = sef.nd_ODE_Static(GetGS).num_v
    ([get_src_edges(gd, i) for i in 1:num_v],
     [get_dst_edges(gd, i) for i in 1:num_v])
end

end #module
