coefvar(xs) = std(xs) / mean(xs)

snr(xs) = mean(xs) / std(xs)

# Significance testing

pvalue(null, x) = (sum(null .> x) + 1) / (length(null) + 1)

"""
    NullDistribution{T<:Real} <: AbstractVector{T}
    NullDistribution(f::Function, ŷ::AbstractArray, y::AbstractArray, n::Int=10_000)

A type representing the null distribution of a test statistic of type `T`. `NullDistribution`s can be constructed by applying function `f` to `ŷ` predictions and `y` ground truths in `n` random permutations of `y`.
"""
struct NullDistribution{T<:Real} <: AbstractVector{T}
    dbn::Vector{T}
end
function NullDistribution(f::Function, ŷ::AbstractArray, y::AbstractArray, n::Int=10_000)
    y = deepcopy(y)
    res = Vector{typeof(f(ŷ,y))}(undef, n)

    for i = 1:n
        shuffle!(y)
        res[i] = f(ŷ, y)
    end

    NullDistribution(res)
end

Base.size(x::NullDistribution) = size(x.dbn)
Base.getindex(x::NullDistribution, i::Integer) = x.dbn[i]
Base.IndexStyle(::NullDistribution) = IndexLinear()
(H0::NullDistribution{T})(x::T) where T<:Real = pvalue(H0, x)
