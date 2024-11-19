using NetworkDynamics
####
#### Diffusion system
####
Base.@propagate_inbounds function diffusionedge!(e, v_s, v_d, _, _)
    e[1] = v_s[1] - v_d[1]
    nothing
end
diffusion_edge() = EdgeModel(; g=AntiSymmetric(diffusionedge!), outdim=1, pdim=0)

Base.@propagate_inbounds function diffusion_dedge!(de, e, v_s, v_d, _, _)
    de[1] = 100.0 * (sin(v_s[1] - v_d[1]) - e[1])
    de[2] = 100.0 * (sin(v_d[1] - v_s[1]) - e[2])
    nothing
end
diffusion_dedge() = EdgeModel(; f=diffusion_dedge!, dim=2, pdim=0, g=Fiducial(dst=1:1, src=2:2))

Base.@propagate_inbounds function diffusionvertex!(dv, _, esum, _, _)
    dv[1] = esum[1]
    nothing
end
diffusion_vertex() = VertexModel(; f=diffusionvertex!, dim=1, pdim=0, g=StateMask(1:1))

####
#### inhomogenious kuramoto system
####
Base.@propagate_inbounds function kuramoto_edge!(e, θ_s, θ_d, (K,), t)
    e[1] = K * sin(θ_s[1] - θ_d[1])
end
static_kuramoto_edge() = EdgeModel(; g=AntiSymmetric(kuramoto_edge!), outdim=1, pdim=1)

Base.@propagate_inbounds function kuramoto_vertex!(dθ, θ, esum, (ω,), t)
    dθ[1] = ω + esum[1]
end
kuramoto_vertex_1d() = VertexModel(; f=kuramoto_vertex!, pdim=1, sym=[:θ], g=StateMask(1:1))

Base.@propagate_inbounds function kuramoto_inertia!(dv, v, esum, (P,), t)
    dv[1] = v[2]
    dv[2] = P - 1.0 * v[2]
    dv[2] += esum[1]
end
kuramoto_vertex_2d() = VertexModel(; f=kuramoto_inertia!, dim=2, pdim=1, sym=[:θ, :ω], g=StateMask(1:1));

####
#### powergrid model
####
Base.@propagate_inbounds function piline!(out_src, out_dst, src, dst, p, t)
    L, R, C1, C2 = p
    ω = 2π*50
    Vsrc = Complex(src[1], src[2])
    Vdst = Complex(dst[1], dst[2])
    imain = (Vsrc - Vdst)/(R + im*ω*L)
    iC1 = Vsrc*im*ω*C1
    iC2 = Vdst*im*ω*C2
    isrc = - (imain + iC1)
    idst = imain - iC2
    out_dst[1] = -real(idst)
    out_dst[2] = -imag(idst)
    out_src[1] = -real(isrc)
    out_src[2] = -imag(isrc)
    nothing
end
piline() = EdgeModel(; g=piline!, outsym=(;src=[:src_i_r, :src_i_i], dst=[:dst_i_r, :dst_i_i]),
                       psym=[:L, :R, :C1, :C2])

Base.@propagate_inbounds function pq!(du, u, isum, (P,Q), t)
    ic = Complex(isum[1], isum[2])
    S = Complex(P, Q)
    uc = S/conj(ic)
    du[1] = real(uc)
    du[2] = imag(uc)
    nothing
end
pqnode() = VertexModel(; f=pq!, sym=[:u_r, :u_i], g=StateMask(1:2), psym=[:P, :Q], mass_matrix=0)

Base.@propagate_inbounds function gen!(dv, v, isum, p, T)
    # unpack parameters
    H, P, D, Ω, E_f, T_d_dash, T_q_dash, X_q_dash, X_d_dash, X_d, X_q = p
    Ω_H = Ω / (2*H)
    # unpack inputs
    i = Complex(isum[1], isum[2])
    # unpack state
    u = Complex(v[1], v[2])
    θ = v[3]
    ω = v[4]

    i_c = 1im*i*exp(-1im*θ)
    e_c = 1im*u*exp(-1im*θ)
    p = real(u * conj(i))
    e_d = real(e_c)
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    dθ = ω
    de_d = (1 / T_q_dash)* (- e_d + (X_q - X_q_dash)* i_q)
    de_q = (1 / T_d_dash)* (- e_q - (X_d - X_d_dash) * i_d + E_f)
    de_c = de_d + 1im*de_q
    du = -1im*de_c*exp(1im*θ)+ u*1im*ω
    dω = (P - D*ω - p- (X_q_dash - X_d_dash)*i_d* i_q)*Ω_H

    # pack du
    dv[1] = real(du)
    dv[2] = imag(du)
    dv[3] = dθ
    dv[4] = dω

    nothing
end
generator() = VertexModel(; f=gen!, sym=[:u_r, :u_i, :θ, :ω], g=StateMask(1:2),
            psym=[:H, :P, :D, :Ω, :E_f, :T_d_dash, :T_q_dash, :X_d_dash, :X_q_dash, :X_d, :X_q])
