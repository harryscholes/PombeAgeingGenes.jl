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
    NullDistribution(dbn::Vector{T}) where T<:Real = new{T}(sort(dbn))
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

abstract type AlternativeHypothesis end
struct Greater <: AlternativeHypothesis end
struct Less <: AlternativeHypothesis end

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

function PValue(instance::T, null::NullDistribution{T};
                alternative::AlternativeHypothesis=Greater()) where T<:Real
    n = length(null)
    moreextreme = if alternative isa Greater
        _pvalue_greater(instance, null, n)
    elseif alternative isa Less
        _pvalue_less(instance, null, n)
    end
    PValue(moreextreme / (n + 1))
end

function PValue(instance::T, null::AbstractVector{T};
                alternative::AlternativeHypothesis=Greater()) where T<:Real
    PValue(instance, NullDistribution(null); alternative=alternative)
end

_pvalue_greater(instance, null, n) = _pvalue_alternative(instance, null, n:-1:1, >)
_pvalue_less(instance, null, n) = _pvalue_alternative(instance, null, 1:n, <)

function _pvalue_alternative(instance, null, iter, op)
    moreextreme = 1

    for i = length(null):-1:1
        if op(instance, null[i])
            break
        end
        moreextreme += 1
    end

    return moreextreme
end
