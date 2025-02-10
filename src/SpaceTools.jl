module SpaceTools

using Dates
using DimensionalData
using DimensionalData.Dimensions
using LinearAlgebra
using Makie
using TimeseriesTools
using Statistics
using Unitful, DimensionfulAngles
using Latexify, UnitfulLatexify
using RollingWindowArrays
using InteractiveViz
using DSP, SignalAnalysis

# export AbstractDataSet, DataSet
export degap, rectify_datetime, resolution, samplingrate, smooth
export timeshift, tnorm
export tplot!, tplot, tplot_panel, tplot_panel!
export tsheat, tlims!, ylabel, plot_attributes, add_labels!
export LMN
export rotate, fac_mat, mva, mva_mat, check_mva_mat
export modify_meta, amap, ω2f

const DD = DimensionalData
const AbstractDimType = Union{AbstractDimStack,AbstractDimArray}
const AbstractDimMatrix = Union{DimensionalData.AbstractDimMatrix,TimeseriesTools.AbstractDimMatrix}
const AbstractDimVector = Union{DimensionalData.AbstractDimVector,TimeseriesTools.AbstractDimVector}

# include("dataset.jl")
include("mhd.jl")
include("methods.jl")
include("timeseries.jl")
include("timeseries/tplot.jl")
include("timeseries/interactive.jl")
include("timeseries/spectrum.jl")
include("utils.jl")
include("utils/dimensiondata.jl")
include("plot.jl")
include("cotrans/coordinate.jl")
include("cotrans/rotate.jl")
include("cotrans/fac.jl")
include("cotrans/mva.jl")

end