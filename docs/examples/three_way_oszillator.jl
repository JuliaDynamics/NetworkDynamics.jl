using WGLMakie
# using GLMakie
# GLMakie.activate!()
using WGLMakie.GeometryBasics
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkit
using ModelingToolkit: D_nounits as D, t_nounits as t
# Function that creates a surface plot with a potential function
function plot_potential_sphere(potential_func::Function; 
                             phi_range=0:π/100:π,      # 0 to π radians
                             theta_range=0:π/100:2π,   # 0 to 2π radians
                             colormap=:viridis)
    # Create the grid of angles
    phis = phi_range
    thetas = theta_range
    
    # Calculate the radius for each point based on the potential function
    radii = [potential_func(phi, theta) for phi in phis, theta in thetas]
    
    # Calculate the 3D coordinates
    # xs = radii .* sin.(phis) .* cos.(thetas')
    # ys = radii .* sin.(phis) .* sin.(thetas')
    # zs = radii .* cos.(phis) .* ones(size(thetas'))
    xs = 1 .* sin.(phis) .* cos.(thetas')
    ys = 1 .* sin.(phis) .* sin.(thetas')
    zs = 1 .* cos.(phis) .* ones(size(thetas'))
    
    # Create the plot
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data)
    surf = surface!(ax, xs, ys, zs; color=radii, colorrange=extrema(radii))
    Colorbar(fig[1, 2], surf)
    
    return fig
end

# Example usage with a constant potential (creates a regular sphere)
plot_potential_sphere((phi, theta) -> 1.0)

# Example with a varying potential
plot_potential_sphere((phi, theta) -> 1.0 + 0.3 * sin(3*theta))

# More complex example with both phi and theta dependence
plot_potential_sphere((phi, theta) -> 1.0 + 0.2 * sin(2*phi) * cos(3*theta))

pot = (N) -> (phi, theta) -> begin
    pot_phi = 1-cos(2*phi)
    pot_theta = 0.5*(1+cos(N*theta))
    0.5*pot_theta*pot_phi
end
plot_potential_sphere(pot(3))


@mtkmodel Oszillator begin
    @structural_parameters begin
        N = 2
    end
    @variables begin
        θ(t)=0 # 0..π
        φ(t)=0 # 0..2π
        ωθ(t)=0 # angular velocity in theta
        ωφ(t)=0 # angular velocity in phi
        pot(t)
        ∂pot∂θ(t)
        ∂pot∂φ(t)
        # φ_force(t)
        # θ_force(t)
    end
    @parameters begin
        M = 1.0
        pot_strength = 0.001
    end
    @equations begin
        pot    ~ pot_strength * (1-cos(2*θ)) * (1 - cos(N*φ))
        ∂pot∂θ ~ pot_strength * (2*sin(2*θ)) * (1 - cos(N*φ))
        ∂pot∂φ ~ pot_strength * (1-cos(2*θ)) * (N*sin(N*φ))
        D(φ) ~ ωφ
        D(θ) ~ ωθ
        M*D(ωφ) ~ -∂pot∂φ
        M*D(ωθ) ~ -∂pot∂θ
    end
end

@mtkbuild oszi = Oszillator(N=2)
prob = ODEProblem(oszi, [oszi.θ => π/4,oszi.φ => π/2, oszi.ωθ => 1.0, oszi.ωφ => 0.0], (0,10))
sol = solve(prob, Rodas5P())

lines(sol, idxs=oszi.θ)
lines(sol, idxs=oszi.φ)

lines(sol, idxs=oszi.ωθ)
lines(sol, idxs=oszi.ωφ)

lines(sol, idxs=oszi.∂pot∂θ)
lines(sol, idxs=oszi.∂pot∂φ)

# Function that creates a surface plot with a potential function
function plot_sol(sol, N)
    pot(θ, φ) = (1-cos(2*θ)) * (1 - cos(N*φ))
    # Create the grid of angles
    thetas = 0:π/100:π
    phis = 0:π/100:2π
    # Calculate the radius for each point based on the potential function
    poti = [pot(theta, phi) for theta in thetas, phi in phis]
    xs = sin.(thetas) .* cos.(phis')
    ys = sin.(thetas) .* sin.(phis')
    zs = cos.(thetas) .* ones(size(phis'))

    # Create the plot
    fig = Figure()
    # ax = Axis3(fig[1, 1], aspect=:data)
    ax = LScene(fig[1,1])
    surf = surface!(ax, xs, ys, zs; color=poti, colorrange=extrema(poti))
    # Colorbar(fig[1, 2], surf)

    pos(θ, φ) = 1.01*Point3f(sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ))
    # Plot the trajectory
    trajectory = [pos(sol(t, idxs=[oszi.θ, oszi.φ])...) for t in range(sol.t[begin], sol.t[end], length=1000)]
    lines!(ax, trajectory, linewidth=2, color=:red)

    return fig
end
plot_sol(sol,2)
