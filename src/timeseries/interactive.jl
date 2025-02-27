using InteractiveViz: Continuous1D, iviz

abstract type AbstractRangeFunction end

### 1D range function that takes a range as input and returns an array
# struct RangeFunction1D{F,L} <: InteractiveViz.Continuous1D
struct RangeFunction1D{F,L} <: AbstractRangeFunction
    f::F
    xmin::L
    xmax::L
end

function sample(data::RangeFunction1D, xrange; samples=10000)
    xmin = first(xrange)
    xmax = last(xrange)
    x, y = data.f((xmin, xmax))

    # Limit to 10000 points if needed
    if length(x) > samples
        @info "Data resampled to $samples points"
        indices = round.(Int, range(1, length(x), length=samples))
        x = x[indices]
        y = ndims(y) == 1 ? y[indices] : y[:, indices]
    end

    (; x, y)
end

limits(data::RangeFunction1D) = (data.xmin, data.xmax, nothing, nothing)

get_xrange(limit) = (limit.origin[1], limit.origin[1] + limit.widths[1])

"""
Remove the resolution-based updates and only update based on axis limit to improve performance
"""
function InteractiveViz.iviz(f, data::RangeFunction1D; delay=DEFAULTS.delay)
    lims = limits(data)
    r = range(lims[1], lims[2]; length=2)
    qdata = sample(data, r)
    x = Observable(qdata.x)
    y = Observable(qdata.y)
    fap = f(x, y)

    if current_axis().limits[] == (nothing, nothing)
        xlims!(current_axis(), lims[1], lims[2])
    end

    ax = current_axis()
    reset_limits!(ax)

    axislimits = ax.finallimits
    prev_xrange = Observable(get_xrange(axislimits[]))

    function update(lims)
        xrange = get_xrange(lims)
        # Update if new range extends beyond previously loaded range
        prev_xmin, prev_xmax = prev_xrange[]
        needs_update = xrange[1] < prev_xmin || xrange[2] > prev_xmax

        # Add range check to avoid unnecessary data fetching
        if needs_update
            qdata = sample(data, xrange)
            x.val = qdata.x
            y[] = qdata.y
            prev_xrange[] = xrange
        end
    end

    on(axislimits) do axlimits
        if @isdefined(redraw_limit)
            close(redraw_limit)
        end
        redraw_limit = Timer(x -> update(axlimits), delay)
    end

    return fap
end

flatten(x) = collect(Iterators.flatten(x))
sample(tas, trange, args...; kwargs...) = flatten(tplot_spec.(tas, trange..., args...; kwargs...))

function iviz_api!(ax::Axis, tas, t0, t1, args...; delay=DEFAULTS.delay, kwargs...)
    specs = Observable(sample(tas, (t0, t1), args...; kwargs...))
    plotlist!(ax, specs)
    reset_limits!(ax)

    axislimits = ax.finallimits
    # Keep track of the previous range
    prev_xrange = Observable(get_xrange(axislimits[]))

    function update(lims)
        xrange = get_xrange(lims)
        # Update if new range extends beyond previously loaded range
        prev_xmin, prev_xmax = prev_xrange[]
        needs_update = xrange[1] < prev_xmin || xrange[2] > prev_xmax

        if needs_update
            trange = x2t.(xrange)
            specs[] = sample(tas, trange, args...; kwargs...)
            prev_xrange[] = xrange
        end
    end

    on(axislimits) do axlimits
        if @isdefined(redraw_limit)
            close(redraw_limit)
        end
        redraw_limit = Timer(x -> update(axlimits), delay)
    end
    return specs[]
end

iviz_api(tas, args...; kwargs...) = iviz_api!(current_axis(), tas, args...; kwargs...)