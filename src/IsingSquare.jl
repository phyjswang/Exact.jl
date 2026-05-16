# 2D Ising model

struct IsingSquare{TJ<:Real} <: AbstractModel
    J::TJ
end

IsingSquare(; J::Real = 1.0) = IsingSquare(J)

"""
Helmholtz free Energy
Infinite-size classical Ising model on the square lattice
H = -J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Exactly solved models in statistical mechanics (2007)]
"""
function getfe(model::IsingSquare, T::Real)
    @unpack J = model
    K = J/T
    twoK = 2K
    sinh2K = sinh(twoK)
    k = inv(sinh2K * sinh2K)
    invk = inv(k)
    cosh2K2 = cosh(twoK)^2
    F(θ) = begin
        s, c = sincos(2θ)
        distance = hypot(muladd(-k, c, 1), k * s)
        -T * (log(2) + log(muladd(invk, distance, cosh2K2))) / (2π)
    end
    return quadgk(F, 0, π)
end
