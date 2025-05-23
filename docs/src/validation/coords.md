# Validate Coordinate Transformations with IRBEM and PySPEDAS

<!-- TODO: gse2gsm -->

See [Coordinate Systems](../explanations/coords.md) for more information.

References: [`cotrans`](@ref), [test_cotrans.py - PySPEDAS](https://github.com/spedas/pyspedas/blob/master/pyspedas/cotrans_tools/tests/cotrans.py)

```@example coords
using PySPEDAS
using SPEDAS
using PythonCall
using DimensionalData
using Chairmarks
```

```@example coords
@py import pyspedas.cotrans_tools.tests.cotrans: CotransTestCases

pytest = CotransTestCases()
pytest.test_cotrans()

tha_pos = PySPEDAS.get_data(DimArray, "tha_pos")
jl_tha_pos = set(tha_pos, Dim{:time}=>Ti)
jl_tha_pos_geo = cotrans(jl_tha_pos', "GEI", "GEO")'
jl_tha_pos_gsm = cotrans(jl_tha_pos', "GEI", "GSM")'
pyspedas.cotrans("tha_pos", "tha_pos_new_gsm", coord_in="GEI", coord_out="GSM")
py_tha_pos_geo = PySPEDAS.get_data(DimArray, "tha_pos_new_geo")
py_tha_pos_gsm = PySPEDAS.get_data(DimArray, "tha_pos_new_gsm")

@assert jl_tha_pos_geo ≈ parent(py_tha_pos_geo)
@assert isapprox(jl_tha_pos_gsm, parent(py_tha_pos_gsm), rtol=1e-3)
```

### Benchmark

Julia's implementation is about 40 times faster than IRBEM's (Fortran) implementation, and 80 times faster than PySPEDAS's (Python) implementation.

```@example coords
@b gei2geo($jl_tha_pos), cotrans($jl_tha_pos', "GEI", "GEO"), pyspedas.cotrans("tha_pos", "tha_pos_new_geo", coord_in="GEI", coord_out="GEO")
```

```@example coords
@b cotrans($jl_tha_pos', "GEI", "GSM"), pyspedas.cotrans("tha_pos", "tha_pos_new_gsm", coord_in="GEI", coord_out="GSM")
```