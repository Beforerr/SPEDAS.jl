@testitem "plot recipes" begin
    using Makie
    @test_nowarn dualplot((rand(3), rand(4)); plotfunc=scatterlines!)
end

@testitem "LinesPlot" begin
    using Makie
    using DimensionalData
    using Unitful
    @test_nowarn linesplot((rand(3), rand(4)))
    ys = [[1, 2, 4] [3, 4, 10]]
    @test_nowarn linesplot(ys)
    @test_nowarn linesplot([10, 20, 30], ys)
    @test_nowarn linesplot([[1, 2, 4], [3, 4, 10, 11]])

    t = Ti(range(DateTime(2000), step=Hour(1), length=4))
    A = rand(t, Y(1:5))
    Au = A * 1u"nT"
    @test_nowarn linesplot(Au)

    @test_throws Unitful.DimensionError let
        f = Figure()
        ax = Axis(f[1, 1]; SpaceTools.axis_attributes(Au)...)
        linesplot!(ax, t.val, Au.data)
        f
    end
end