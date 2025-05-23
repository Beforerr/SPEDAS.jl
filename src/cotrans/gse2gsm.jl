"""
    gse2gsm_mat(time)

Compute the GSE to GSM transformation matrix.
"""
function gse2gsm_mat(time)
    dipole_gei = geo2gei_mat(time) * cdipdir(time)
    # Get sun direction parameters
    _, _, sra, sdec, obliq = csundir(time)
    # Calculate sun direction vector in GEI
    sun_gei = SA[cos(sra)*cos(sdec), sin(sra)*cos(sdec), sin(sdec)]
    # Calculate ecliptic pole direction in GEI
    pole_gei = SA[0.0, -sin(obliq), cos(obliq)]

    # Calculate cross product of dipole and sun directions
    gmgs = cross(dipole_gei, sun_gei)
    # Calculate magnitude of cross product
    rgmgs = norm(gmgs)
    # Calculate cosine and sine of GSE-GSM angle
    cdze = dot(pole_gei, dipole_gei) / rgmgs
    sdze = dot(pole_gei, gmgs) / rgmgs
    return SA[1 0 0; 0 cdze sdze; 0 -sdze cdze]
end

"""
    gse2gsm(pos, time)

Transform position vector from GSE to GSM coordinates.

The GSE to GSM transformation is given by the matrix T3 = <- psi,X>, where the rotation angle psi 
is the GSE-GSM angle. This transformation is a rotation in the GSE YZ plane from the 
GSE Z axis to the GSM Z axis.

See also: [`gse2gsm_mat`](@ref)
"""
gse2gsm(pos, time) = gse2gsm_mat(time) * pos

"""
    gsm2gse(pos, time)

Transform position vector from GSM to GSE coordinates.

See also: [`gse2gsm`](@ref)
"""
gsm2gse(pos, time) = inv(gse2gsm_mat(time)) * pos