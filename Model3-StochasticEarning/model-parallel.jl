#Per-period Utility Function
@everywhere function uc(c::Float64; σ=2.0)
    return (c^(1-σ))/(1-σ)
end
