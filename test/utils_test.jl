using NetworkDynamics

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

@testset "test utils" begin
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
        using NetworkDynamics: compg
        a = Lib.diffusion_edge()
        b = Lib.diffusion_edge()
        @test compg(a) == compg(b)
        a = Lib.diffusion_edge_closure()
        b = Lib.diffusion_edge_closure()
        @test compg(a) != compg(b)

        fid = Lib.diffusion_edge_fid()
        ode = Lib.diffusion_odeedge()

        odevert = Lib.diffusion_vertex()

        kura_edge   = Lib.kuramoto_edge()
        kura_second = Lib.kuramoto_second()
    end

    @testset "find_identical" begin
        using NetworkDynamics: _find_identical

        v1 = Lib.kuramoto_second()
        @test _find_identical(v1, 1:10) == [collect(1:10)]
        v2 = Lib.diffusion_vertex()
        v3 = Lib.kuramoto_first()

        # v2 and v3 are equal when it comes to the function!!
        vs = [v1,v2,v3,v2,v2,v1,v1,v3]

        @test _find_identical(vs, eachindex(vs)) == [[1,6,7],[2,4,5],[3,8]]

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

    @testset "unique_mappings" begin
         using NetworkDynamics: unique_mappings
         a = [1,2,3]
         b = [1,2,3]
         @test unique_mappings(a, b) == Dict(1=>1, 2=>2, 3=>3)
         a = [:a, :b, :c]
         b = [1,2,3]
         @test unique_mappings(a, b) == Dict(:a=>1, :b=>2, :c=>3)
         a = [:a, :b, :c]
         b = [1,1,3]
         @test unique_mappings(a, b) == Dict(:a=>1, :b=>1, :c=>3)
         a = [:a, :a, :c]
         b = [1,1,3]
         @test unique_mappings(a, b) == Dict(:c=>3)
         a = [:a, :a, :c]
         b = [1,1,-3]
         @test unique_mappings(abs, a, b) == Dict(:c=>3)
    end

    @testset "flatrange to compbine ranges" begin
        using NetworkDynamics: flatrange
        r1 = 1:10
        r2 = 11:12
        @test flatrange((src=r1,dst=r2)) == 1:12
        @test flatrange(r1) == r1
        r3 = 12:13
        @test_throws AssertionError flatrange((src=r1,dst=r3))
    end

    @testset "significant digits print" begin
        using NetworkDynamics: str_significant
        @test str_significant(0.0; sigdigits=3) == "0"
        @test str_significant(1.23456789; sigdigits=3) == "1.23"
        @test str_significant(9.81; sigdigits=3) == "9.81"
        @test str_significant(9.81; sigdigits=2) == "9.8"
        @test str_significant(1002.1; sigdigits=3) == "1e3"
        @test str_significant(0.00000122; sigdigits=3) == "1.22e-6"
        @test str_significant(-123.294191; sigdigits=5) == "-123.29"
    end

    @testset "batch stride" begin
        using NetworkDynamics: BatchStride, _fullstride, _fullrange, _range
        using Static

        b1 = BatchStride(2, 3)
        b2 = BatchStride(2, (3, 1))

        @test _fullstride(b1) == 3
        @test _fullstride(b2) == 4
        Test.@inferred _fullstride(b1)
        Test.@inferred _fullstride(b2)

        @test _fullrange(b1, 2) == 2:(2*3+2-1)
        @test _fullrange(b2, 2) == 2:(2*(3+1)+2-1)

        @test _range(b1, 1) == 2:4
        @test _range(b1, 2) == 5:7
        @test _range(b2, 1, 1) == 2:4
        @test _range(b2, 1, 2) == 5:5
        @test _range(b2, 2, 1) == 6:8
        @test _range(b2, 2, 2) == 9:9

        # named batch stride
        bn = BatchStride(2, (;src=3, dst=1))
        @test _fullstride(bn) == 4
        Test.@inferred _fullstride(bn)
        @test _fullrange(bn, 2) == 2:(2*(3+1)+2-1)
        @test _range(bn, 1, 1) == 2:4
        @test _range(bn, 1, 2) == 5:5
        @test _range(bn, 2, 1) == 6:8
        @test _range(bn, 2, 2) == 9:9
    end
end
