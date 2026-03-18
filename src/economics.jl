# Helpful economic functions used in many models 

## Utility
@inline _crra(c::Real, γ::Real) = γ == 1.0 ? _logutility(c) : _crra_strict(c, γ)
@inline _powerutility(c::Real, γ::Real) = γ == 1.0 ? _logutility(c) : _powerutility_strict(c, γ)
@inline _crra_prime(c::Real, γ::Real) = c^(-γ)
@inline _logutility(c::Real) = log(c)
@inline _logutility_prime(c::Real) = 1/c
@inline _crra_strict(c::Real, γ::Real) = (c^(1-γ) - 1)/(1-γ)
@inline _powerutility_strict(c::Real, γ::Real) = c^(1-γ)/(1-γ)

"CRRA utility function with risk aversion parameter `γ`."
const crra = (u = _crra, u′ = _crra_prime)

"Log utility function."
const logutility = (u = _logutility, u′ = _logutility_prime)

"CRRA utility function with risk aversion parameter `γ ≠ 1`."
const crra_strict = (u = _crra_strict, u′ = _crra_prime)

"Power utility function (CRRA without `-1` adjustment) with risk aversion parameter `γ`."
const powerutility = (u = _powerutility, u′ = _crra_prime)

"Power utility function (CRRA without `-1` adjustment) with risk aversion parameter `γ`."
const powerutility_strict = (u = _powerutility_strict, u′ = _crra_prime)

## CES
"CES aggregator over vector `x` with weights vector `w` and elasticity of substitution `ε`."
function ces(x, w, ε::Real)
    if ε == 1.0
        return cobb_douglas(x, w)
    else
        return ces_strict(x, w, ε)
    end
end

function ces_strict(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, ε::Real) 
    a = (ε - 1) / ε
    
    # Initialize accumulator with the correct type to maintain type stability
    s = zero(promote_type(eltype(x), eltype(w))) 
    
    @inbounds @simd for i in eachindex(w, x)
        s += w[i] * x[i]^a
    end
    
    return s^(1/a)
end

@inline function ces_strict(x::SVector{<:Real}, w::SVector{<:Real}, ε::Real) 
    a = (ε - 1) / ε
    return sum(w .* x.^a)^(1/a)
end

"CES aggregator over two inputs `x₁`, `x₂` with weight `w` on `x₁` and `1-w` on `x₂` and elasticity of substitution `ε`."
@inline function ces(x₁::Real, x₂::Real, w::Real, ε::Real)
    ε == 1.0 ? cobb_douglas(x₁, x₂, w) : ces_strict(x₁, x₂, w, ε)
end
@inline function ces_strict(x₁::Real, x₂::Real, w::Real, ε::Real)
    a = (ε - 1) / ε
    return (w * x₁^a + (1 - w) * x₂^a)^(1 / a)
end

"Marginal product of input `i` in a CES function with inputs `x`, weights `w`, and elasticity of substitution `ε`."
@inline function ces_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, ε::Real, i::Integer)
    Y = ces(x, w, ε)
    return ces_prime(x, w, ε, Y, i)
end

"Marginal product of input `i` in a CES function with inputs `x`, weights `w`, elasticity of substitution `ε`, and aggregate Y given."
@inline function ces_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, ε::Real, Y::Real, i::Integer)
    return ces_prime(x[i], w[i], ε, Y) 
end

"Marginal product in a CES function of input `xᵢ` with weight `wᵢ`, elasticity of substitution `ε`, and aggregate Y given."
@inline function ces_prime(xi::Real, w::Real, ε::Real, Y::Real)
    return w * (Y / xi)^(1/ε) 
end
        

## Cobb-Douglas
"Cobb-Douglas aggregator over vector `x` with exponents vector `w`."
function cobb_douglas(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}) 
    # Initialize accumulator with the correct type to maintain type stability
    p = one(promote_type(eltype(x), eltype(w))) 
    
    @inbounds @simd for i in eachindex(w, x)
        p *= x[i]^w[i]
    end
    
    return p
end

@inline function cobb_douglas(x::SVector{<:Real}, w::SVector{<:Real}) 
    return prod(x.^w)
end

"Cobb-Douglas aggregator over two inputs `x₁`, `x₂` with exponents `w` and `1-w`."
@inline cobb_douglas(x₁::Real, x₂::Real, w::Real) = x₁^w * x₂^(1 - w)

"Marginal product of input `i` in a Cobb-Douglas function with inputs `x` and exponents `w`."
@inline function cobb_douglas_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, i::Integer)
    Y = cobb_douglas(x, w)
    return cobb_douglas_prime(x, w, Y, i)
end

"Marginal product of input `i` in a Cobb-Douglas function with inputs `x` and exponents `w` with aggregate Y given."
@inline function cobb_douglas_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, Y::Real, i::Integer)
    return cobb_douglas_prime(x[i], w[i], Y)
end

"Marginal product in a Cobb-Douglas function of inputs `xᵢ` and exponent `wᵢ` with aggregate Y given."
@inline function cobb_douglas_prime(xi::Real, wi::Real, Y::Real)
    return wi * Y / xi
end

"Marginal product in a Cobb-Douglas function of inputs `xᵢ, x₂` and exponent `wᵢ`."
@inline function cobb_douglas_prime(x₁::Real, x₂::Real, w::Real, Y::Real, i::Integer)
    if i == 1
        return w * Y / x₁
    elseif i == 2
        return (1-w) * Y / x₂
    else
        throw(ArgumentError("Input index i must be 1 or 2 for two-input Cobb-Douglas. Got i=$i"))
    end
end

"Marginal product in a Cobb-Douglas function of inputs `xᵢ, x₂` and exponent `wᵢ` with aggregate Y given."
@inline function cobb_douglas_prime(x₁::Real, x₂::Real, w::Real,  i::Integer)
    Y = cobb_douglas(x₁, x₂, w)
    return cobb_douglas_prime(x₁, x₂, w, Y, i)
end

"Marginal products of inputs in a Cobb-Douglas function with inputs `x` and exponents `w` with aggregate Y given."
@inline function cobb_douglas_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}, Y::Real)
    return cobb_douglas_prime.(x, w, Y)
end

"Marginal products of inputs in a Cobb-Douglas function with inputs `x` and exponents `w`."
@inline function cobb_douglas_prime(x::AbstractArray{<:Real}, w::AbstractArray{<:Real})
    Y = cobb_douglas(x, w)
    return cobb_douglas_prime(x, w, Y::Real)
end