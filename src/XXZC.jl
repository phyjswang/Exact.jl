# XXZ Chain
# H = J∑ᵢ(SˣᵢSˣᵢ₊₁ + SʸᵢSʸᵢ₊₁) + ΔJ∑ᵢSᶻᵢSᶻᵢ₊₁

struct XXZC <: AbstractModel
    J::Real
    Δ::Real
end

"""
Ground state energy
Infinite-size XXZ chain
H = J∑ᵢ(SˣᵢSˣᵢ₊₁ + SʸᵢSʸᵢ₊₁) + ΔJ∑ᵢSᶻᵢSᶻᵢ₊₁

[Ref: V. E. Korepin, N. M. Bogoliubov, A. G. Izergin, Quantum inverse scattering method and correlation functions, Cambridge University Press (1997)]
"""
function getgse(model::XXZC)
    @unpack J, Δ = model
    η = 1/2*acos(-Δ)
    rslt, err = quadgk(λ -> J*sin(2η)^2 / (π - 2η) /(cosh(π*λ/(π-2η)) * cosh(λ+1.0im*η) * cosh(λ - 1.0im * η)), -Inf, Inf)
    return J*Δ/4 - rslt, err
end
