# Transverse Field Ising Chain
# H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

struct TFIC{TJ<:Real,Tg<:Real} <: AbstractModel
    J::TJ
    g::Tg
end

@inline function _tfic_energy(J::Real, g::Real, k::Real)
    s, c = sincos(k)
    halfJ = J / 2
    return hypot(muladd(halfJ, c, g), halfJ * s) / 2
end

"""
Ground state energy
Infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getgse(model::TFIC)
    @unpack J, g = model
    halfJ = J / 2
    scale = g + halfJ
    parameter = 2 * J * g / (scale * scale)
    return -abs(scale) * ellipe(parameter) / π, 0
end

"""
Helmholtz free energy
Infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getfe(model::TFIC, T::Real)
    @unpack J, g = model
    return quadgk(k -> -T / π * _log2cosh(_tfic_energy(J, g, k) / T), 0, π)
end

"""
Internal energy
Infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getie(model::TFIC, T::Real)
    @unpack J, g = model
    return quadgk(k -> begin
        y = _tfic_energy(J, g, k)
        -y * tanh(y / T) / π
    end, 0, π)
end

"""
Loschmidt rate function
Infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: PRB 87, 195104 (2013)]
"""
function getlrf(m0::TFIC, m1::TFIC, t::Real)
    J = m0.J
    J ≠ m1.J && error("J must be the same for both models by our convention.")
    g0 = m0.g
    g1 = m1.g
    g0 *= 2/J
    g1 *= 2/J
    J /= 4

    ε(g::Real, s::Real, c::Real) = 2J * hypot(g - c, s)

    function θ(g::Real, s::Real, c::Real)
        θ1 = atan(s, g - c) / 2
        if θ1 < 0
            return θ1 + π/2
        else
            return θ1
        end
    end

    y = quadgk(k -> begin
        s, c = sincos(k)
        φ = θ(g0, s, c) - θ(g1, s, c)
        sinφ, cosφ = sincos(φ)
        weight = muladd(sinφ * sinφ, cis(-2 * t * ε(g1, s, c)), cosφ * cosφ)
        -log(weight) / π
    end, 0, π)

    return real(y[1]), y[2]
end
