using DimensionalData.Lookups
using HybridArrays

const timeDimType = (DimensionalData.TimeDim, Dim{:time})

abstract type FrequencyDim{T} <: Dimension{T} end
@dim 𝑓 FrequencyDim "Frequency"

function hybridify(A, dims)
    sizes = ntuple(ndims(A)) do i
        i in dims ? StaticArrays.Dynamic() : size(A, i)
    end
    HybridArray{Tuple{sizes...}}(A)
end

function hybridify(A::AbstractDimArray, dim)
    rebuild(A, hybridify(parent(A), dim))
end

hybridify(A; query=nothing) = 
    hybridify(A, dimnum(A, something(query, TimeDim)))

function standardize(x::AbstractDimArray; floatify=true)
    # Convert integer values to floats
    floatify && eltype(x) <: Integer && (x = modify(float, x))
    # Check if any of the dimensions match our time dimension types
    x = any(d -> d isa Dim{:time}, dims(x)) ? set(x, Dim{:time} => Ti) : x
end
tdim(t) = Ti(t)
tdim(t::DD.Dimension) = t

function tvec(A::AbstractDimArray; query=nothing)
    dim = timedim(A, query)
    DimArray(vec(parent(A)), dim; metadata=meta(A), name=name(A))
end

"""
    nt2ds(nt_arr, dim; fields=propertynames(first(nt_arr)))

Convert a NamedTuple array to a DimStack of DimArrays.
"""
function nt2ds(nt_arr, dim; fields=propertynames(first(nt_arr)))
    DimStack([
        DimArray(getfield.(nt_arr, field), dim; name=field)
        for field in fields
    ])
end

function nt2ds(nt_arr; sym=:time)
    dim = Dim{sym}(getfield.(nt_arr, sym))
    # filter the time dimension
    fields = propertynames(first(nt_arr))
    fields = filter(field -> field != sym, fields)
    nt2ds(nt_arr, dim; fields)
end

function rename(da::AbstractDimArray, new_name)
    rebuild(da; name=new_name)
end

function _dimnum(x, query = nothing)
    # If query is nothing, choose the time dimension or the first dimension
    if isnothing(query)
        hasdim(x, TimeDim) ? dimnum(x, TimeDim) : 1
    else
        dimnum(x, query)
    end
end
    
"""
    amap(f, a, b)

Apply a function `f` to the intersection of `a` and `b`.

https://github.com/rafaqz/DimensionalData.jl/issues/914
"""
function amap(f, a::AbstractDimArray, b::AbstractDimArray)
    shared_selectors = DimSelectors(a)[DimSelectors(b)]
    f(a[shared_selectors], b[shared_selectors])
end