# https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.cotrans_tools.cotrans_lib.subgeo2gei

@inline function gei2geo_mat(time)
    gst = calculate_gst_alt(time)
    cgst = cos(gst)
    sgst = sin(gst)
    return SA[cgst sgst 0; -sgst cgst 0; 0 0 1]
end

geo2gei_mat(time) = inv(gei2geo_mat(time))

"""
    gei2geo(pos, time)

Converts geocentric equatorial inertial (GEI) coordinates to geographical (GEO) coordinates.
"""
gei2geo(pos, time) = gei2geo_mat(time) * pos

"""
    geo2gei(pos, time)

Converts geographical (GEO) coordinates to geocentric equatorial inertial (GEI) coordinates.
"""
geo2gei(pos, time) = geo2gei_mat(time) * pos