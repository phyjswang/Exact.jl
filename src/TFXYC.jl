# TFXY Chain
# H = -∑ᵢ(JₓSˣᵢSˣᵢ₊₁ + JySʸᵢSʸᵢ₊₁) - h∑ᵢSᶻᵢ

struct TFXYC <: AbstractModel
    Jx::Real
    Jy::Real
    h::Real
end

"""
Ground state energy
infinite-size TFXY chain
H = -∑ᵢ(JxSˣᵢSˣᵢ₊₁ + JySʸᵢSʸᵢ₊₁) - h∑ᵢSᶻᵢ

[Ref: SciPost Phys. Lect. Notes 82 (2024)]
"""
function getgse(model::TFXYC)
    @unpack Jx, Jy, h = model
    Jx /= 4
    Jy /= 4
    h /= 2
    J = Jx+Jy
    κ = (Jx-Jy)/(Jx+Jy)
    ε(k) = 2*abs(J)*√((h/J - cos(k))^2 + κ^2 * sin(k)^2)
    return quadgk(k -> -ε(k)/(2π), 0, π)
end
