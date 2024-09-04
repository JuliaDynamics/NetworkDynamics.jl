using NetworkDynamics

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "test utils" begin
    @testset begin
        using InteractiveUtils
        for t in subtypes(NetworkDynamics.Coupling)
            println(NetworkDynamics._styled_coupling(t()))
        end
    end

    @testset "subscript" begin
        using NetworkDynamics: subscript
        @test subscript(10) == "₁₀"
        @test subscript(5) == "₅"
    end

    @testset "style symbol array" begin
        using NetworkDynamics: stylesymbolarray
        syms = [:a, :b, :c]
        defaults = [1, 2, nothing]
        stylesymbolarray(syms, defaults, Dict(1 => :red, 2 => :orange))
        stylesymbolarray(syms, defaults, Dict(1 => :red, 2 => :red))
        NetworkDynamics.ND_FACES
    end

    @testset "Test Component Library" begin
        using NetworkDynamics: compf
        a = Lib.diffusion_edge()
        b = Lib.diffusion_edge()
        @test compf(a) == compf(b)
        a = Lib.diffusion_edge_closure()
        b = Lib.diffusion_edge_closure()
        @test compf(a) != compf(b)

        fid = Lib.diffusion_edge_fid()
        ode = Lib.diffusion_odeedge()

        odevert = Lib.diffusion_vertex()
        odevert_const = Lib.diffusion_vertex_constraint()

        kura_edge   = Lib.kuramoto_edge()
        kura_second = Lib.kuramoto_second()
    end

    @testset "find_identical" begin
        using NetworkDynamics: _find_identical

        v1 = Lib.kuramoto_second()
        @test _find_identical(v1, 1:10) == [collect(1:10)]
        v2 = Lib.diffusion_vertex()
        v3 = Lib.diffusion_vertex_constraint()

        # v2 and v3 are equal when it comes to the function!!
        vs = [v1,v2,v3,v2,v2,v1,v1,v3]

        @test _find_identical(vs, eachindex(vs)) == [[1,6,7],[2,3,4,5,8]]

        es = [Lib.diffusion_edge(),
              Lib.diffusion_edge_closure(),
              Lib.diffusion_edge_closure(),
              Lib.diffusion_edge_fid()]
        @test _find_identical(es, eachindex(es)) == [[1], [2], [3], [4]]
    end

    @testset "algin strings" begin
        using NetworkDynamics: align_strings

        align_strings(["1 &=> 2", ":foobar &=> 8"])

        align_strings(["1 &=> 2 & k->4",
                       ":foobar &=> 8"])

        align_strings(["1 &=> 2 & k &-> 4",
                       ":foobar &=> 8 & :a &-> 4"])

        align_strings([":a &=> 2 & k &-> 4",
                       ":foobar &=> 8 & :a &-> 4"])

        align_strings(["row &with\n&line break", "second &row"])
    end
end
