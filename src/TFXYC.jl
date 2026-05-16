# TFXY Chain
# H = -∑ᵢ(JₓSˣᵢSˣᵢ₊₁ + JySʸᵢSʸᵢ₊₁) - h∑ᵢSᶻᵢ

struct TFXYC{TJx<:Real,TJy<:Real,Th<:Real} <: AbstractModel
    Jx::TJx
    Jy::TJy
    h::Th
end

"""
Ground state energy
Infinite-size TFXY chain
H = -∑ᵢ(JxSˣᵢSˣᵢ₊₁ + JySʸᵢSʸᵢ₊₁) - h∑ᵢSᶻᵢ

[Ref: SciPost Phys. Lect. Notes 82 (2024)]
"""
function getgse(model::TFXYC)
    @unpack Jx, Jy, h = model
    Jx /= 4
    Jy /= 4
    h /= 2
    J = Jx + Jy
    κ = (Jx - Jy) / J
    absJ = abs(J)
    invJ = inv(J)
    ε(k) = begin
        s, c = sincos(k)
        2 * absJ * hypot(h * invJ - c, κ * s)
    end
    return quadgk(k -> -ε(k) / (2π), 0, π)
end
