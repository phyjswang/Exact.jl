module Exact

using UnPack
using SpecialFunctions
using QuadGK

abstract type AbstractModel end

export TFIC
export getgse
export getfe
export getie

include("TFIC.jl")

end
