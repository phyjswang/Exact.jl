# Transverse Field Ising Chain
# H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

struct TFIC <: AbstractModel
    J::Real
    g::Real
end

"""
Ground state energy
infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getgse(model::TFIC)
    @unpack J, g = model
    return -2/π * (g+J/2) * ellipe(2J/g / (1 + J/2g)^2) / 2
end

"""
Helmholtz free energy
infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getfe(model::TFIC, T::Real)
    @unpack J, g = model
    y = quadgk(k -> -T * 1/π * log(cosh(1/2/T*g*√(1+(J/2g)^2 + J/g * cos(k)))), 0, π)
    return -T * log(2) + y[1], y[2]
end

"""
Internal energy
infinite-size transverse field Ising chain
H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - g Σᵢ Sˣᵢ

[Ref: Annals of Physics 57, 79 (1970)]
"""
function getie(model::TFIC, T::Real)
    @unpack J, g = model
    y(k) = 1/2*g*√(1+(J/2g)^2 + J/g * cos(k))
    return quadgk(k ->T^2 * 1/π *  -y(k) / T^2 * tanh(y(k)/T), 0, π)
end

"""
Loschmidt rate function
infinite-size transverse field Ising chain
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

    ε(g::Real,k::Real) = 2J* √((g-cos(k))^2 + sin(k)^2)

    function θ(g::Real,k::Real)
        θ1 = atan(sin(k) / (g-cos(k)))/2
        if θ1 < 0
            return θ1 + π/2
        else
            return θ1
        end
    end

    φ(g0::Real,g1::Real,k::Real) = θ(g0,k) - θ(g1,k)

    y = quadgk(k -> - 2 / 2π * log(cos(φ(g0,g1,k))^2 + sin(φ(g0,g1,k))^2 * exp(-2*1.0im*t*ε(g1,k))), 0, π)

    return real(y[1]), y[2]
end
