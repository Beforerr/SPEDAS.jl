# Coordinate Systems

This package defines common coordinate systems used in heliophysics and space physics research.

```@contents
Pages = ["coords.md"]
Depth = 2:3
```

## Coordinate Transformations

```@docs; canonical=false
rotate
```

A comprehensive description of the transformations can be found in [hapgoodSpacePhysicsCoordinate1992](@citet)

## Standard Coordinate Systems

* GSE (Geocentric Solar Ecliptic)
* GSM (Geocentric Solar Magnetic)

Systems based on the Earth's rotation axis

* GEO (Geographic)
* GEI (Geocentric Equatorial Inertial)
* J2000

Systems based on the dipole axis of the Earth's magnetic field

* SM (Solar Magnetic)
* MAG (Geomagnetic)

More information can be found in the the following links

- [https://stereo-ssc.nascom.nasa.gov/coordinates_explanation.shtml](https://stereo-ssc.nascom.nasa.gov/coordinates_explanation.shtml)
- [Geocentric systems](https://www.mssl.ucl.ac.uk/grid/iau/extra/local_copy/SP_coords/geo_sys.htm]

### Coordinate transformations between geocentric systems

```@docs; canonical=false
cotrans
```

References:
- [Coordinate transformations between geocentric systems](https://www.mssl.ucl.ac.uk/grid/iau/extra/local_copy/SP_coords/geo_tran.htm)

`SPEDAS.jl` implements native Julia functions

1. [`geo2gei`](@ref) and [`gei2geo`](@ref): Geographic (GEO) to Geocentric Equatorial Inertial (GEI) coordinate transformations and vice versa. 
2. [`gse2gsm`](@ref) and [`gsm2gse`](@ref): Geocentric Solar Ecliptic (GSE) to Geocentric Solar Magnetic (GSM) coordinate transformations and vice versa.

For other coordinate systems, including magnetic coordinate calculations, we leverage [IRBEM-LIB](https://prbem.github.io/IRBEM/) through its Julia interface [IRBEM.jl](https://github.com/Beforerr/IRBEM.jl).

```@example coords
using DimensionalData
using Speasy, SPEDAS
using CairoMakie

pos_gse = get_data("cda/THC_L1_STATE/thc_pos_gse", "2015-10-16", "2015-10-17") |> DimArray

pos_gsm = cotrans(pos_gse', "GSM")'
pos_sm = cotrans(pos_gse', "SM")'
pos_geo = cotrans(pos_gse', "GEO")'

tplot((pos_gse, pos_gsm, pos_sm, pos_geo))
```

## Specialized Coordinate Systems

The package also provides transformations for analysis-specific coordinate systems:

### Field-Aligned Coordinates (FAC)

A local coordinate system defined relative to the ambient magnetic field direction, useful for studying plasma waves and particle distributions.

```@docs; canonical=false
fac_mat
```

### Minimum Variance Analysis (MVA) and Boundary Normal Coordinates (LMN)

A coordinate system derived from the eigenvalues and eigenvectors of the magnetic field variance matrix, commonly used in analyzing current sheets, discontinuities, and wave propagation.

References:

- [Minimum and Maximum Variance Analysis](https://ui.adsabs.harvard.edu/abs/1998ISSIR...1..185S)
- [https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.minvar](https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.minvar)

```@docs; canonical=false
mva_eigen
mva
check_mva_eigen
```

Error estimates for MVA:

```@docs; canonical=false
SPEDAS.Δφij
SPEDAS.B_x3_error
```

See also: [Comparison with PySPEDAS](../validation/pyspedas.md#minimum-variance-analysis).

## Reference

- [The International Radiation Belt Environment Modeling (IRBEM) library](https://prbem.github.io/IRBEM/)


- [PySPEDAS](https://pyspedas.readthedocs.io/en/latest/coords.html)
- [geopack](https://github.com/tsssss/geopack): Python version of geopack and Tsyganenko models
- [geospacelab](https://github.com/JouleCai/geospacelab): A python-based library to collect, manage, and visualize geospace data (e.g. OMNI, geomagnetic indices, EISCAT, DMSP, SWARM, TEC, AMPERE, etc.).
- [SPEDAS Cotrans](https://spedas.org/wiki/index.php?title=Cotrans)