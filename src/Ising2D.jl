# 2D Ising model

struct Ising2D <: AbstractModel
    J::Real
end

Ising2D(;J::Real = 1.0) = Ising2D(J)

"""
Helmholtz free Energy
Infinite-size classical Ising chain on the square lattice
H = -J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Exactly solved models in statistical mechanics (2007)]
"""
function getfe(model::Ising2D, T::Real)
    @unpack J = model
    K = J/T
    k = 1/(sinh(2*K))^2
    F(θ) = -log(2*(cosh(2*K)^2+k^-1 *√(1+k^2-2*k*cos(2*θ)))) / (2*π) * T
    return quadgk(F, 0, π)
end
