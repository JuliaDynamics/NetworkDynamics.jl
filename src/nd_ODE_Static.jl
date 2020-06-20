module nd_ODE_Static_mod

using ..NetworkStructures
using ..NDFunctions
using ..Utilities


export nd_ODE_Static
export StaticEdgeFunction

#=  =#

# If the type of dx and x and the e_buffer match, we swap out v_array for x
@inline function prep_gls(dx::T, x::T, gls, e_array::T) where {T}
    gls[1].graph_data.gdb.v_array = x
    # Does not allocate, all graph layers are assumed to refer to the same
    # buffer gdb.
    gls
end

# If the type of dx and x do not match, we swap initialize a new GraphData object
# that is based on the type of dx for the edge buffer.
# Some solvers take the derivative with respect to time, thus x will not be dual
# but dx will be, leading to errors otherwise
@inline function prep_gls(dx, x, gls, e_array)
    new_e_array = similar(dx, length(e_array))
    (GraphLayer(gl.edges!,
        gl.aggregator,
        gl.graph,
        gl.graph_structure,
        GraphData(x, new_e_array, gl.graph_structure)) for gl in gls)
end



@Base.kwdef struct nd_ODE_Static{T1, GL}
    vertices!::T1
    graph_layers::GL
    parallel::Bool # enables multithreading for the core loop
end

function (d::nd_ODE_Static)(dx, x, p, t)

    gls = prep_gls(dx, x, d.graph_layers, d.graph_layers[1].graph_data.gdb.e_array)

    p_offset = 0

    for gl in gls
        gs = gl.graph_structure
        es! = gl.edges!
        gd = gl.graph_data

        @nd_threads d.parallel for i in (1:gs.num_e) .+ p_offset
            es!.f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i), t)
        end

        p_offset += gs.num_e
    end

    let gs = gls[1].graph_structure,  gd = gls[1].graph_data

        @nd_threads d.parallel for i in 1:gs.num_v
            maybe_idx(d.vertices!,i).f!(view(dx,gs.v_idx[i]), gd.v[i], p_v_idx(p, i), t, Tuple{Float64}(gl.aggregator(gl.graph_data.e_s_v[i], gl.graph_data.e_d_v[i]) for gl in  gls)...)
        end
    end

    nothing
end


function (d::nd_ODE_Static)(x, p, t, ::Type{GetGD})

    gls = prep_gls(dx, x, d.graph_layers, d.graph_layers[1].graph_data.gdb.e_array)

    p_offset = 0

    for gl in gls
        gs = gl.graph_structure
        es! = gl.edges!
        gd = gl.graph_data

        @nd_threads d.parallel for i in 1:gs.num_e
            maybe_idx(es!, i).f!(gd.e[i], gd.v_s_e[i], gd.v_d_e[i], p_e_idx(p, i + p_offset), t)
        end

        p_offset += gs.num_e
    end

    println("Warning, only returning first layer data. API Change needed.")
    iterate(gls)[1].graph_data
end


function (d::nd_ODE_Static)(::Type{GetGS})
    d.graph_layers[1].graph_structure
end


# For compatibility with PowerDynamics

struct StaticEdgeFunction
    nd_ODE_Static
end

function (sef::StaticEdgeFunction)(x, p, t)
    gd = sef.nd_ODE_Static(x, p, t, GetGD)

    (gd.e_s_v, gd.e_d_v)
end

end #module
