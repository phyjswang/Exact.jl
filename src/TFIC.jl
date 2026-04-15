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
    return -T * (log(2) + 1/π * quadgk(k -> log(cosh(1/2/T*g*√(1+(J/2g)^2 + J/g * cos(k)))), 0, π)[1])
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
    return T^2 * 1/π * quadgk(k -> -y(k) / T^2 * tanh(y(k)/T), 0, π)[1]
end
