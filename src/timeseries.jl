using RollingWindowArrays

function tnorm(da)
    norm.(eachrow(da))
end

function timeshift(ta; dim=1, t0=nothing)
    td = dims(ta, dim)
    times = td.val.data
    t0 = something(t0, times[1])

    new_dim_name = Symbol("Time after ", t0)
    new_dim = Dim{new_dim_name}(times .- t0)

    DimArray(ta.data, (new_dim, otherdims(ta, dim)...), name=ta.name, metadata=ta.metadata)
end

function resolution(times; tol=2)
    dt = diff(times)
    dt0 = eltype(dt)(1)
    dtf_mean, relerr = mean_relerr(dt ./ dt0)
    if relerr > exp10(-tol - 1)
        @warn "Time resolution is is not approximately constant (relerr ≈ $relerr)"
    end
    round(Integer, dtf_mean) * dt0
end

resolution(da::AbstractDimType; dim=Ti, kwargs...) =
    resolution(dims(da, dim).val; kwargs...)

samplingrate(da) = NoUnits(1u"s" / resolution(da))

function smooth(da::AbstractDimArray, span::Integer; dims=Ti, suffix="_smoothed", kwargs...)
    new_da = mapslices(da; dims) do slice
        mean.(RollingWindowArrays.rolling(slice, span; kwargs...))
    end
    rebuild(new_da; name=Symbol(da.name, suffix))
end

function degap(da::DimArray; dim=Ti)
    dims = otherdims(da, dim)
    rows = filter(x -> !any(isnan, x), eachslice(da; dims))
    if !isempty(rows)
        cat(rows...; dims)
    else
        similar(da, (0, size(da, 2)))
    end
end

function rectify_datetime(da; tol=2, kwargs...)
    times = dims(da, Ti)
    t0 = times[1]
    dtime = Quantity.(times.val .- t0)
    new_times = TimeseriesTools.rectify(Ti(dtime); tol)[1]
    set(da, Ti => new_times .+ t0)
end