"""
    calculate_gst(time)

Calculate Greenwich sidereal time (GST) in radians from the given time.

Reference: https://aa.usno.navy.mil/faq/GAST, https://github.com/JuliaAstro/AstroLib.jl/blob/main/src/ct2lst.jl
"""
function calculate_gst(time)
    # ct2lst returns local sidereal time in hours (0-24)
    # For Greenwich, longitude=0, so GST = LST
    return ct2lst(0, jdcnv(time)) * 2π / 24
end

"""
    calculate_gst_alt(time)

Alternative implementation of Greenwich sidereal time calculation based on the
reference algorithm from `pyspedas.cotrans_tools.csundir_vect`.
"""
function calculate_gst_alt(time::DateTime)
    # Extract time components
    iyear = year(time)
    idoy = dayofyear(time)
    # Julian day and Greenwich mean sidereal time calculation
    fday = Time(time).instant / Day(1)
    jj = 365 * (iyear - 1900) + floor((iyear - 1901) / 4) + idoy
    dj = jj - 0.5 + fday
    gst = mod(279.690983 + 0.9856473354 * dj + 360.0 * fday + 180.0, 360.0)
    return gst * 2π / 360
end

"""
    csundir(time)

Calculate the direction of the sun.

# Returns
- `gst`: Greenwich mean sidereal time (radians)
- `slong`: Longitude along ecliptic (radians)
- `sra`: Right ascension (radians)
- `sdec`: Declination of the sun (radians)
- `obliq`: Inclination of Earth's axis (radians)

# Notes
- This function calculates various parameters related to the sun's position
  needed for coordinate transformations.
"""
function csundir(time)
    # Convert time to year, day of year, hour, minute, second
    dt = time isa DateTime ? time : DateTime(time)
    year = Dates.year(dt)
    doy = Dates.dayofyear(dt)
    hour = Dates.hour(dt)
    minute = Dates.minute(dt)
    second = Dates.second(dt)

    # Julian day and Greenwich mean sidereal time
    pisd = π / 180.0
    fday = (hour * 3600.0 + minute * 60.0 + second) / 86400.0
    jj = 365 * (year - 1900) + floor((year - 1901) / 4) + doy
    dj = jj - 0.5 + fday
    gst = mod(279.690983 + 0.9856473354 * dj + 360.0 * fday + 180.0, 360.0) * pisd

    # Longitude along ecliptic
    vl = mod(279.696678 + 0.9856473354 * dj, 360.0)
    t = dj / 36525.0
    g = mod(358.475845 + 0.985600267 * dj, 360.0) * pisd
    slong = (vl + (1.91946 - 0.004789 * t) * sin(g) + 0.020094 * sin(2.0 * g)) * pisd

    # Inclination of Earth's axis
    obliq = (23.45229 - 0.0130125 * t) * pisd
    sob = sin(obliq)
    cob = cos(obliq)

    # Aberration due to Earth's motion around the sun (about 0.0056 deg)
    pre = (0.005686 - 0.025e-4 * t) * pisd

    # Declination of the sun
    slp = slong - pre
    sind = sob * sin(slp)
    cosd = sqrt(1.0 - sind^2)
    sc = sind / cosd
    sdec = atan(sc)

    # Right ascension of the sun
    sra = π - atan((cob / sob) * sc, -cos(slp) / cosd)

    return gst, slong, sra, sdec, obliq
end
