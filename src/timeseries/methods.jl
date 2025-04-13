function tmean(x; dims=timedim(x))
    mean(x; dims)
end

function tmedian(x; dims=timedim(x))
    nanmedian(x; dims)
end

"""
    tsubtract(x, f=nanmedian; dims=timedim(x))

Subtract a statistic (default function `f`: `nanmedian`) along dimensions (default: time dimension) from `x`.
"""
function tsubtract(x, f=nanmedian; dims=timedim(x))
    x .- f(parent(x); dims=dimnum(x, dims))
end

function tnorm(x; dims=Ti)
    norm.(eachslice(x; dims))
end

"""
    proj(a, b)

Vector projection of a vector `a` on (or onto) a nonzero vector `b`.

# References: [Wikipedia](https://en.wikipedia.org/wiki/Vector_projection)

See also: [`sproj`](@ref), [`oproj`](@ref)
"""
function proj(a, b)
    b̂ = normalize(ustrip(b))
    return dot(a, b̂) * b̂
end

"""Scalar projection"""
function sproj(a, b)
    return dot(a, b) / norm(b)
end

function tsproj(a, b; dims=Ti)
    sproj.(eachslice(a; dims), eachslice(b; dims))
end

function tproj(a, b; dims=Ti)
    proj.(eachslice(a; dims), eachslice(b; dims))
end

"""Vector rejection"""
oproj(a, b) = a - proj(a, b)

function toproj(v, B; dims=Ti)
    res = oproj.(eachslice(v; dims), eachslice(B; dims))
    tstack(res)
end

"""

References:
- https://docs.xarray.dev/en/stable/generated/xarray.cross.html
"""
function tcross(x, y; dims=Ti, stack=nothing)
    stack = @something stack (ndims(x) == 2)
    res = cross.(eachslice(x; dims), eachslice(y; dims))
    stack ? tstack(res) : res
end

function tdot(x, y; dims=Ti)
    dot.(eachslice(x; dims), eachslice(y; dims))
end

norm_combine(x; dims=1) = cat(x, norm.(eachslice(x; dims)); dims=setdiff(1:ndims(x), dims))

"""
    tnorm_combine(x; dims=timedim(x), name=:magnitude)

Calculate the norm of each slice along `query` dimension and combine it with the original components.
"""
function tnorm_combine(x; dims=timedim(x), name=:magnitude)
    data = norm_combine(parent(x); dims=dimnum(x, dims))

    # Replace the original dimension with our new one that includes the magnitude
    odim = otherdims(x, dims) |> only
    odimType = basetypeof(odim)
    new_odim = odimType(vcat(odim.val, name))
    new_dims = map(d -> d isa odimType ? new_odim : d, DD.dims(x))
    return rebuild(x, data, new_dims)
end


for f in (:tsubtract,)
    @eval $f(args...; kwargs...) = x -> $f(x, args...; kwargs...)
end