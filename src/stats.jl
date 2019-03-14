coefvar(xs) = std(xs) / mean(xs)

snr(xs) = mean(xs) / std(xs)

# Significance testing

"""
    NullDistribution([f::Function,] ŷ::AbstractArray, y::AbstractArray[, n=10_000])

Construct a null distribution for probabilities `ŷ` against ground truths `y`.

`f` is applied to `n` random permutations of `y`. If `f` is not supplied PR AUC is calculated.

# Example

```julia
ŷ = rand(1000)
y = rand(Bool, 1000)
nd = NullDistribution(ŷ, y)
pvalue(nd, auc(PR(ŷ, y)))
```
"""
struct NullDistribution{T<:Real} <: AbstractVector{T}
    dbn::Vector{T}
end

function NullDistribution(f::Function, ŷ::AbstractVector, y::AbstractVector, n=10_000)
    y = deepcopy(y)
    dbn = map(_->f(ŷ, shuffle!(y)), 1:n)
    NullDistribution(dbn)
end

function NullDistribution(ŷ::AbstractVector, y::AbstractVector, n=10_000)
    f(ŷ, y) = auc(PR(ŷ, y))
    NullDistribution(f, ŷ, y, n)
end

Base.size(x::NullDistribution) = size(x.dbn)
Base.getindex(x::NullDistribution, i::Integer) = x.dbn[i]
Base.IndexStyle(::NullDistribution) = IndexLinear()

# p-values

function _pvalue(dbn::AbstractVector{<:Real}, instance::Real)
    (sum(dbn .> instance) + 1) / (length(dbn) + 1)
end

"""
    PValue([f::Function,] ŷ::AbstractArray, y::AbstractArray[, n=10_000])

Calucluate a p-value for the probabilities `ŷ` against ground truths `y`.

`f` is applied to `n` random permutations of `y`. If `f` is not supplied PR AUC is calculated.

# Example

```julia
ŷ = rand(1000)
y = rand(Bool, 1000)
PValue(ŷ, y)
```
"""
struct PValue{T<:Real}
    p::T
end

function PValue(f::Function, ŷ::AbstractVector, y::AbstractVector, n=10_000)
    PValue(NullDistribution(f, ŷ, y, n), f(ŷ, y))
end

function PValue(ŷ::AbstractVector, y::AbstractVector, n=10_000)
    f(ŷ, y) = auc(PR(ŷ, y))
    dbn = NullDistribution(f, ŷ, y, n)
    instance = f(ŷ, y)
    PValue(_pvalue(dbn, instance))
end

PValue(dbn::NullDistribution, instance) = PValue(_pvalue(dbn, instance))
