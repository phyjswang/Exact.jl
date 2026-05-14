# 2D Ising model

struct IsingSquare <: AbstractModel
    J::Real
end

IsingSquare(;J::Real = 1.0) = IsingSquare(J)

"""
Helmholtz free Energy
Infinite-size classical Ising model on the square lattice
H = -J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Exactly solved models in statistical mechanics (2007)]
"""
function getfe(model::IsingSquare, T::Real)
    @unpack J = model
    K = J/T
    k = 1/(sinh(2*K))^2
    F(θ) = -log(2*(cosh(2*K)^2+k^-1 *√(1+k^2-2*k*cos(2*θ)))) / (2*π) * T
    return quadgk(F, 0, π)
end
