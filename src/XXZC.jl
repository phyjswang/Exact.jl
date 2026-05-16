# XXZ Chain
# H = J∑ᵢ(SˣᵢSˣᵢ₊₁ + SʸᵢSʸᵢ₊₁) + ΔJ∑ᵢSᶻᵢSᶻᵢ₊₁

struct XXZC{TJ<:Real,TΔ<:Real} <: AbstractModel
    J::TJ
    Δ::TΔ
end

"""
Ground state energy
Infinite-size XXZ chain
H = J∑ᵢ(SˣᵢSˣᵢ₊₁ + SʸᵢSʸᵢ₊₁) + ΔJ∑ᵢSᶻᵢSᶻᵢ₊₁

[Ref: V. E. Korepin, N. M. Bogoliubov, A. G. Izergin, Quantum inverse scattering method and correlation functions, Cambridge University Press (1997)]
"""
function getgse(model::XXZC)
    @unpack J, Δ = model
    η = acos(-Δ) / 2
    width = π - 2η
    prefactor = J * sin(2η)^2 / width
    cosη2 = cos(η)^2
    rslt, err = quadgk(λ -> begin
        sinhλ = sinh(λ)
        denominator = cosh(π * λ / width) * muladd(sinhλ, sinhλ, cosη2)
        prefactor / denominator
    end, -Inf, Inf)
    return J*Δ/4 - rslt, err
end
