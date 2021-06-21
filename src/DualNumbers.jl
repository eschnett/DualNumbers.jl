module DualNumbers

"""
Dual numbers with a scalar dual part
"""
struct Dual{T} <: Real
    primal::T
    dual::T
    Dual{T}(x, y) where {T} = new{T}(T(x), T(y))
end
Dual(x::T1, y::T2) where {T1,T2} = Dual{promote_type(T1, T2)}(x, y)
Dual{T}(x) where {T} = Dual{T}(x, 0)
Dual(x) = Dual(x, 0)
export Dual

# Only compare primal parts
Base.:(==)(x::Dual, y::Dual) = x.primal == y.primal
Base.:(<)(x::Dual, y::Dual) = x.primal < y.primal

Base.:(==)(x::Dual, y::Real) = x == Dual(y)
Base.:(<)(x::Dual, y::Real) = x < Dual(y)
Base.:(==)(x::Real, y::Dual) = Dual(x) == y
Base.:(<)(x::Real, y::Dual) = Dual(x) < y

Base.promote_rule(::Type{Dual{T1}}, ::Type{Dual{T2}}) where {T1,T2} = Dual{promote_type(T1, T2)}
Base.promote_rule(::Type{T1}, ::Type{Dual{T2}}) where {T1<:Real,T2} = Dual{promote_type(T1, T2)}

Base.zero(x::Dual) = Dual(zero(x.primal))
Base.one(x::Dual) = Dual(one(x.primal))

Base.:+(x::Dual) = Dual(+x.primal, +x.dual)
Base.:-(x::Dual) = Dual(-x.primal, -x.dual)

Base.:+(x::Dual, y::Dual) = Dual(x.primal + y.primal, x.dual + y.dual)
Base.:-(x::Dual, y::Dual) = Dual(x.primal - y.primal, x.dual - y.dual)

Base.inv(x::Dual) = Dual(inv(x.primal), -x.dual / x.primal^2)

Base.:*(x::Dual, y::Dual) = Dual(x.primal * y.primal, x.dual * y.primal + x.primal * y.dual)
Base.:/(x::Dual, y::Dual) = Dual(x.primal / y.primal, x.dual / y.primal - x.primal / y.primal^2 * y.dual)
Base.:\(x::Dual, y::Dual) = Dual(x.primal \ y.primal, x.primal \ y.dual - (x.primal^2 \ y.primal) * x.dual)

Base.literal_pow(::typeof(^), x::Dual, ::Val{0}) = one(x)
function Base.literal_pow(::typeof(^), x::Dual, ::Val{n}) where {n}
    return Dual(Base.literal_pow(^, x.primal, Val(n)), n * Base.literal_pow(^, x.primal, Val(n - 1)) * x.dual)
end
Base.:^(x::Dual, n::Integer) = n == 0 ? one(x) : Dual(x.primal^n, n * x.primal^(n - 1) * x.dual)
Base.:^(x::Dual, y::Real) = (@assert !(y isa Dual); Dual(x.primal^y, y * x.primal^(y - 1) * x.dual))
Base.:^(x::Dual, y::Dual) = error("not implemented")

Base.cos(x::Dual) = Dual(cos(x.primal), -x.dual * sin(x.primal))
Base.log(x::Dual) = Dual(log(x.primal), x.dual / x.primal)
Base.sin(x::Dual) = Dual(sin(x.primal), x.dual * cos(x.primal))
Base.sqrt(x::Dual) = Dual(sqrt(x.primal), x.dual / (2 * sqrt(x.primal)))

################################################################################

"""
    derivative(f, x)
    derivative(f, xs, i)

Evaluate the partial derivative of the scalar-argument function `f(x)`
at `x`, or of the array-argument function `f(xs)` at `xs[i]`.
"""
derivative(f, x::T) where {T<:Real} = f(Dual(x, 1)).dual
function derivative(f, xs::AbstractArray{T}, n::Integer) where {T<:Real}
    xs′ = similar(xs, Dual{T})
    for i in 1:length(xs′)
        xs′[i] = Dual(xs[i], i == n)
    end
    return f(xs′).dual
end
export derivative

end
