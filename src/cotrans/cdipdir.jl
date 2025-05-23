"""
    cdipdir(time)

Compute dipole direction in GEO coordinates.

# Notes
- Computes geodipole axis direction from International Geomagnetic Reference
  Field (IGRF-14) model for time interval 1965 to 2030.
- For time out of interval, computation is made for nearest boundary.

References
- https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.cotrans_tools.cotrans_lib.cdipdir
"""
function cdipdir(time)
    # Convert time to year and day of year
    dt = time isa DateTime ? time : DateTime(time)
    year = Dates.year(dt)
    doy = Dates.dayofyear(dt)

    # Get IGRF-14 parameters
    minyear, maxyear, ga, ha, dg, dh = set_igrf_params()

    # Find base year for IGRF coefficients
    y = year - (year % 5)
    # Ensure y is within the valid range
    y = max(min(y, maxyear - 5), minyear)
    # Do not interpolate beyond boundaries
    year = max(min(year, maxyear), minyear)

    year0 = y
    year1 = y + 5
    maxind = maxyear - 5

    if year1 <= maxind
        f2 = (year + (doy - 1) / 365.25 - year0) / 5.0
        f1 = 1.0 - f2
        g = @. ga[year0] * f1 + ga[year1] * f2
        h = @. ha[year0] * f1 + ha[year1] * f2
    else
        f3 = year + (doy - 1) / 365.25 - maxind
        g = @. ga[year0] + dg * f3
        h = @. ha[year0] + dh * f3
    end


    # Schmidt normalization for spherical harmonic coefficients
    # IGRF model uses Schmidt semi-normalized spherical harmonic coefficients
    s = 1.0
    for i in 2:14
        mn = floor(Int, i * (i - 1) / 2 + 1)
        s = floor(Int, s * (2 * i - 3) / (i - 1))
        g[mn+1] *= s
        h[mn+1] *= s
        g[mn] *= s
        h[mn] *= s
        p = s
        for j in 2:i-1
            aa = 1.0
            if j == 2
                aa = 2.0
            end
            p = p * sqrt(aa * (i - j + 1) / (i + j - 2))
            mnn = mn + j - 1
            g[mnn+1] *= p
            h[mnn+1] *= p
            g[mnn] *= p
            h[mnn] *= p
        end
    end

    g10 = -g[2]  # Adjusting for 1-based indexing in Julia
    g11 = g[3]
    h11 = h[3]

    sq = g11^2 + h11^2
    sqq = sqrt(sq)
    sqr = sqrt(g10^2 + sq)
    s10 = -h11 / sqq
    c10 = -g11 / sqq
    st0 = sqq / sqr
    ct0 = g10 / sqr

    return SA[st0*c10, st0*s10, ct0]
end