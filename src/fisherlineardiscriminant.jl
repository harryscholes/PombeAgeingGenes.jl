#=
Fisher's linear discriminant analysis
=#
using Statistics, LinearAlgebra, StatsBase, MultivariateStats

function probability(M::Discriminant, X::AbstractMatrix)
    y = evaluate(M, X)
    y .= 1 ./ (exp.(y .* -1) .+ 1)
end

# FisherLinearDiscriminant

"""
    FisherLinearDiscriminant(Xp, Xn)

Ref: Pattern Recognition and Machine Learning - Bishop (2006) Section 4.1.4
http://read.pudn.com/downloads190/ebook/893343/Kernel_Methods_for_Pattern_Analysis/0521813972.pdf
"""
struct FisherLinearDiscriminant{T<:Real} <: Discriminant{T}
    w::Vector{T}
    μ::Vector{T}
end
function FisherLinearDiscriminant(Xp::AbstractMatrix, Xn::AbstractMatrix)
    μp = mean(Xp, dims=2)                     # feature means
    μn = mean(Xn, dims=2)
    Zp = Xp .- μp                             # mean centre
    Zn = Xn .- μn
    Cp = Zp * Zp'                             # covariance matrix
    Cn = Zn * Zn'
    Sw = Cp + Cn                              # within-class covariance matrix
    w = cholesky(Sw) \ (μp - μn)              # projection weight vector
    μ = mean(hcat(Xp, Xn), dims=2)            # mean of all data
    FisherLinearDiscriminant(vec(w), vec(μ))
end

function StatsBase.fit(::Type{FisherLinearDiscriminant}, Xp, Xn)
    FisherLinearDiscriminant(Xp, Xn)
end

function kthfold(T::Type{FisherLinearDiscriminant}, X, y, tX; kwargs...)
    M = fit(T, X[:, y .== 1], X[:, y .== 0])
    probability(M, tX)
end

function MultivariateStats.evaluate(M::FisherLinearDiscriminant, X::AbstractVecOrMat)
    y = (X .- M.μ)' * M.w
end
