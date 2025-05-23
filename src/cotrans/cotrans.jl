import IRBEM
using AstroLib: ct2lst, jdcnv

export cotrans

include("coordinate.jl")
include("igrf.jl")
include("csundir.jl")
include("cdipdir.jl")
include("geo2gei.jl")
include("gse2gsm.jl")

"""
    cotrans(A, in, out)
    cotrans(A, out; in=get_coord(A))

Transform the data to the `out` coordinate system from the `in` coordinate system.

This function is a wrapper for `IRBEM.transform`.
"""
function cotrans(A, in, out)
    time = parent(times(A))
    data = IRBEM.transform(time, parent(A), in, out)
    set_coord(rebuild(A; data=data), out)
end

cotrans(A, out; in=get_coord(A)) = cotrans(A, in, out)
cotrans(A, f::Function; dims=1) = map(f, eachslice(parent(A); dims), times(A))

for f in (:gei2geo, :geo2gei, :gse2gsm, :gsm2gse)
    @eval @inline function $f(A)
        dims = dimnum(A, TimeDim)
        data = stack($f, eachslice(parent(A); dims), times(A); dims)
        rebuild(A, data)
    end
    @eval export $f
end
