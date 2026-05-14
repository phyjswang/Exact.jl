module Exact

using UnPack
using SpecialFunctions
using QuadGK
using HCubature

abstract type AbstractModel end

export TFIC
export XXZC
export TFXYC
export IsingSquare
export IsingKagome

export getgse
export getfe
export getie
export getlrf

include("TFIC.jl")
include("XXZC.jl")
include("TFXYC.jl")
include("IsingSquare.jl")
include("IsingKagome.jl")

end
