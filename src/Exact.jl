module Exact

using UnPack
using SpecialFunctions
using QuadGK
using HCubature

abstract type AbstractModel end

@inline function _log2cosh(x::Real)
    ax = abs(x)
    return ax + log1p(exp(-2 * ax))
end

@inline function _cosh3_minus_sinh3(x::Real, c::Real, s::Real)
    return exp(-x) * (c * c + c * s + s * s)
end

@inline _sinh_cosh_difference(x::Real, c::Real, s::Real) = c * s * exp(-x)

@inline _ising_lattice_sum(x, y) = cos(x) + cos(y) - cos(x - y)

export TFIC
export XXZC
export TFXYC
export IsingSquare
export IsingKagome
export IsingTriangular

export getgse
export getfe
export getie
export getlrf

include("TFIC.jl")
include("XXZC.jl")
include("TFXYC.jl")
include("IsingSquare.jl")
include("IsingKagome.jl")
include("IsingTriangular.jl")

end
