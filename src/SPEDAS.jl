module SPEDAS

using Dates
using NanoDates
using DimensionalData
using DimensionalData.Dimensions
using DimensionalData: AbstractDimMatrix, AbstractDimVector
using LinearAlgebra
import NaNMath
using Makie
using Intervals: Interval
using StructArrays
using StaticArrays
using OhMyThreads
using Statistics
using Unitful, DimensionfulAngles
using Latexify, UnitfulLatexify
using RollingWindowArrays
using InteractiveViz
using FFTW, DSP, SignalAnalysis
using Tullio
using Bumper
using InverseFunctions
using NamedTupleTools

export AbstractProduct, Product, SpeasyProduct
export dropna, rectify_datetime, resolution, samplingrate, smooth, tsplit
export timerange, TimeRange
export tclip, timeshift, tmean, tnorm, norm_combine, tnorm_combine, tmean, tcross, tdot
export proj, tproj, sproj, tsproj, oproj, toproj
export tstack, tinterp, tinterp_nans, resample, tresample, tfilter
export find_spikes, replace_outliers
export fill_gaps
export tplot!, tplot, tplot_panel, tplot_panel!
export tsheat, tlims!, tlines!, add_labels!
export ylabel, plot_attributes
export LMN
export rotate, select_rotate, fac_mat, tfac_mat, mva, mva_mat, check_mva_mat
export get_coord, get_coords, set_coord
export standardize, modify_meta, amap, ω2f
export reciprocal_vector, reciprocal_vectors, lingradest
export volumetric_tensor, tetrahedron_quality
export ConstantVelocityApproach, CVA, ConstantThicknessApproach, CTA, DiscontinuityAnalyzer, DA
export Elsässer, σ_c
export spectral_matrix, wavpol, twavpol, wpol_helicity, polarization

const DD = DimensionalData
const AbstractDimType = Union{AbstractDimStack,AbstractDimArray}
const MatrixLike = Union{AbstractArray{<:AbstractVector},AbstractMatrix}

include("DataModel/DataModel.jl")
include("projects/project.jl")
include("mhd.jl")
include("methods.jl")
include("timeseries.jl")
include("timeseries/interactive.jl")
include("timeseries/spectrum.jl")
include("timeseries/methods.jl")
include("timeseries/spike.jl")
include("timeseries/gap.jl")
include("utils.jl")
include("utils/timerange.jl")
include("utils/dimensiondata.jl")
include("types.jl")
include("resampling/resample.jl")
include("resampling/interp.jl")
include("meta.jl")
include("plot/types.jl")
include("plot/transform.jl")
include("plot/attributes.jl")
include("plot/tplot.jl")
include("plot/methods.jl")
include("cotrans/coordinate.jl")
include("cotrans/rotate.jl")
include("cotrans/fac.jl")
include("cotrans/mva.jl")
include("multispacecraft/reciprocal_vector.jl")
include("multispacecraft/tetrahedron.jl")
include("multispacecraft/lingradest.jl")
include("multispacecraft/timing.jl")
include("waves/polarization.jl")
include("waves/helicty.jl")
include("waves/spectral_matrix.jl")

end