# 2D Ising model

struct IsingKagome <: AbstractModel
    J::Real
end

IsingKagome(;J::Real = 1.0) = IsingKagome(J)

"""
Helmholtz free Energy
Infinite-size classical Ising model on the kagome lattice
H = -J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Entropy 27, 799 (2010)]
"""
function getfe(model::IsingKagome, T::Real)
    @unpack J = model
    c = cosh(2J/T)
    s = sinh(2J/T)
    F(x) = -T * log(16*((c^3-s^3)^2+3*c^2)+32*s*(s*c^2-c*s^2) * (cos(x[1])-cos(x[1]-x[2])+cos(x[2]))) *3 / (24*π^2)
    return hcubature(F, [0,0], [2π,2π])
end
