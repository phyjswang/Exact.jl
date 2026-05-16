# 2D Ising model

struct IsingTriangular{TJ<:Real} <: AbstractModel
    J::TJ
end

IsingTriangular(; J::Real = 1.0) = IsingTriangular(J)

"""
Helmholtz free Energy
Infinite-size classical Ising model on the triangular lattice
H = J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Entropy 27, 799 (2010)]
"""
function getfe(model::IsingTriangular, T::Real)
    @unpack J = model
    x = 2J / T
    c = cosh(x)
    s = sinh(x)
    c3ms3 = _cosh3_minus_sinh3(x, c, s)
    F(k) = -T * log(c3ms3 - s * _ising_lattice_sum(k[1], k[2])) / (8π^2)
    rslt = hcubature(F, [0,0], [2π,2π])
    return -T * log(2) + rslt[1], rslt[2]
end
