module Exact

using UnPack
using SpecialFunctions
using QuadGK

abstract type AbstractModel end

export TFIC
export XXZC
export TFXYC
export Ising2D

export getgse
export getfe
export getie
export getlrf

include("TFIC.jl")
include("XXZC.jl")
include("TFXYC.jl")
include("Ising2D.jl")

end
