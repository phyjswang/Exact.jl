# 2D Ising model

struct IsingKagome{TJ<:Real} <: AbstractModel
    J::TJ
end

IsingKagome(; J::Real = 1.0) = IsingKagome(J)

"""
Helmholtz free Energy
Infinite-size classical Ising model on the kagome lattice
H = J ∑<ᵢ,ⱼ> SᵢSⱼ

[Ref: Entropy 27, 799 (2010)]
"""
function getfe(model::IsingKagome, T::Real)
    @unpack J = model
    x = 2J / T
    c = cosh(x)
    s = sinh(x)
    c2 = c * c
    c3ms3 = _cosh3_minus_sinh3(x, c, s)
    sc2mcs2 = _sinh_cosh_difference(x, c, s)
    offset = 16 * (c3ms3 * c3ms3 + 3 * c2)
    amplitude = 32 * s * sc2mcs2
    F(k) = -T * log(muladd(amplitude, _ising_lattice_sum(k[1], k[2]), offset)) / (8π^2)
    return hcubature(F, [0,0], [2π,2π])
end
