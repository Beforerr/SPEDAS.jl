# https://github.com/MakieOrg/Makie.jl/blob/master/src/basic_recipes/series.jl
# https://github.com/rafaqz/DimensionalData.jl/blob/main/ext/DimensionalDataMakie.jl
import Makie: convert_arguments, plot!

struct NoDimConversion <: Makie.ConversionTrait end

@recipe LinesPlot begin
    labels = nothing
    # Makie.MakieCore.documented_attributes(Lines)...
    # resample = 10000
end

Makie.conversion_trait(::Type{<:LinesPlot}) = NoDimConversion()

function Makie.convert_arguments(::Type{<:LinesPlot}, x::AbstractVector, ys::AbstractMatrix)
    curves = map(i -> (x, view(ys, :, i)), 1:size(ys, 2))
    return (curves,)
end

function Makie.convert_arguments(T::Type{<:LinesPlot}, ys::AbstractMatrix)
    Makie.convert_arguments(T, 1:size(ys, 1), ys)
end

"""Convert the vector into a single-column matrix"""
function Makie.convert_arguments(T::Type{<:LinesPlot}, ys::AbstractVector{<:Number})
    return Makie.convert_arguments(T, reshape(ys, :, 1))
end

"""Convert the vector of vectors into a single vector of curves"""
function Makie.convert_arguments(T::Type{<:LinesPlot}, ys::Union{Tuple,AbstractVector})
    tuples = Makie.convert_arguments.(T, ys)
    curves_vec = first.(tuples)
    curves = reduce(vcat, curves_vec)
    return (curves,)
end

function Makie.convert_arguments(::Type{<:LinesPlot}, da::DimensionalData.AbstractDimMatrix)
    x = lookup(dims(da, 1)).data
    labels = SpaceTools.labels(da)
    map(eachcol(parent(da)), labels) do y, label
        S.Lines(x, y; label)
    end
end

function Makie.convert_arguments(::Type{<:LinesPlot}, da::DimensionalData.AbstractDimVector)
    x = lookup(dims(da, 1)).data
    S.Lines(x, parent(da))
end

function Makie.plot!(plot::LinesPlot)
    curves = plot[1]
    nseries = length(curves[])
    for i in 1:nseries
        positions = lift(c -> c[i], plot, curves)
        x = lift(x -> x[1], positions)
        y = lift(x -> x[2], positions)
        lines!(plot, x, y)
    end
end

"""
    linesplot(gp, ta)

Plot a multivariate time series on a panel
"""
function linesplot(gp::Drawable, ta; axis=(;), add_title=DEFAULTS.add_title, kwargs...)
    ax = Axis(gp; axis_attributes(ta; add_title)..., axis...)
    plots = linesplot!(ax, ta; kwargs...)
    PanelAxesPlots(gp, AxisPlots(ax, plots))
end

# function linesplot!(ax, xs, vs::Observable; labels, kwargs...)
#     nseries = size(vs[], 2)
#     plots = []
#     for i in 1:nseries
#         y_col = @lift($vs[:, i])
#         plot = lines!(ax, xs, y_col; label=labels[i], kwargs...)
#         push!(plots, plot)
#     end
#     plots
# end