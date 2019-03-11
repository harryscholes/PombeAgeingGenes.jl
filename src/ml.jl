# CV

"""
    @crossvalidate(X, y[, k=5[, f=stratifiedkfolds]], ex)

Estimate the performance of a model using features `X` and labels `y` in `k` folds using
partitioning function `f` and loop body `ex`.

The variables `Xtrain`, `ytrain`, `Xtest` and `ytest` are available to the cross-validation
Expr `ex`. `ex` should return the predicted probabilities for `Xtest`. `Performance`
prediction metrics for each fold  and a `PR` precision-recall curve are returned.
"""
macro crossvalidate(args...)
    # Required
    X = args[1]
    y = args[2]
    ex = args[end]

    # Optional
    k = 5
    f = stratifiedkfolds

    na = length(args)

    if na == 4
        k = args[3]
    elseif na == 5
        k = args[3]
        f = args[4]
    else
        na != 3 && throw(ArgumentError("wrong number of arguments to @crossvalidate"))
    end

    esc(quote
        local data = shuffleobs(($X, $y))
        local ys = []
        local ŷs = []

        for ((Xtrain, ytrain), (Xtest, ytest)) = ($f)(data, $k)
            ŷ = $ex
            push!(ys, ytest)
            push!(ŷs, ŷ)
        end

        Performance.(ŷs, ys), PR(vcat(ŷs...), vcat(ys...))
    end)
end

# Stratified k-fold cross-validation

function stratifiedkfolds(f, data, k::Integer=5, obsdim=default_obsdim(data))
    data_shuf = shuffleobs(data, obsdim=obsdim)
    lm = MLDataUtils.labelmap(eachtarget(f, data_shuf, obsdim=obsdim))
    val_indices = map(_->Int[], 1:k)

    for (class, indices) = lm
        # Validation indicies corresponding to the values in the labelmap dict
        _, lm_val_indices = kfolds(length(indices), k)
        for i = 1:k
            # Map lm_val_indices indicies back to data indicies
            append!(val_indices[i], lm[class][lm_val_indices[i]])
        end
    end

    train_indices = map(t->setdiff(1:nobs(data; obsdim=obsdim), t), val_indices)
    FoldsView(data_shuf, train_indices, val_indices, obsdim)
end

function stratifiedkfolds(data, k::Integer=5, obsdim=default_obsdim(data))
    stratifiedkfolds(identity, data, k, obsdim)
end

# Performance metrics

# methods for tp, fn, fp, tn
for (f, gt, pr) = ((:tp,true,true), (:fn,true,false), (:fp,false,true), (:tn,false,false))
    f_ = string(f)
    @eval begin
        """
            $($f_)(ŷ::AbstractArray{Bool}, y::AbstractArray{Bool})

        $($f_) for `ŷ` predicted labels against `y` ground truth labels.
        """
        @inline function ($f)(ŷ::AbstractArray{Bool}, y::AbstractArray{Bool})
            ŷ = ŷ[:]; y = y[:]
            count(x->x == $pr, ŷ[y .== $gt])
        end
    end
end

"""
    @validate ex

If `ex` evaluates to `NaN`, return zero, else return the evaluation.
"""
macro validate(ex)
    esc(quote
        local x = $ex
        isnan(x) ? zero(x) : x
    end)
end

import Base: precision
import HypothesisTests: FisherExactTest

"""
    accuracy(TP::Integer, FN::Integer, FP::Integer, TN::Integer)
    accuracy(ŷ::AbstractArray, y::AbstractArray)
    accuracy(c::ConfusionMatrix)

Accuracy = (TP + TN) / (TP + TN + FP + FN)
"""
function accuracy(TP::Integer, FN::Integer, FP::Integer, TN::Integer)
    @validate (TP + TN) / (TP + TN + FP + FN)
end

"""
    precision(TP::Integer, FP::Integer)
    precision(ŷ::AbstractArray, y::AbstractArray)
    precision(c::ConfusionMatrix)

Precision = TP / (TP + FP)
"""
function precision(TP::Integer, FP::Integer)
    @validate TP / (TP + FP)
end

"""
    recall(TP::Integer, FN::Integer)
    recall(ŷ::AbstractArray, y::AbstractArray)
    recall(c::ConfusionMatrix)

Recall = TP / (TP + FN)
"""
function recall(TP::Integer, FN::Integer)
    @validate TP / (TP + FN)
end

"""
    f1(TP::Integer, FN::Integer, FP::Integer)
    f1(ŷ::AbstractArray, y::AbstractArray)
    f1(c::ConfusionMatrix)

F1 = 2 * (Precision * Recall) / (Precision + Recall)
"""
function f1(TP::Integer, FN::Integer, FP::Integer)
    Precision, Recall = precision(TP, FP), recall(TP, FN)
    @validate 2 * (Precision * Recall) / (Precision + Recall)
end

"""
    mcc(TP::Integer, FN::Integer, FP::Integer, TN::Integer)
    mcc(ŷ::AbstractArray, y::AbstractArray)
    mcc(c::ConfusionMatrix)

Matthews correlation coefficient

MCC = ((TP*TN) - (FP*FN)) / √((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
"""
function mcc(TP::Integer, FN::Integer, FP::Integer, TN::Integer)
    @validate ((TP*TN) - (FP*FN)) / √((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
end

# confusion matrix

"""
    ConfusionMatrix

Confusion matrix with shape:

        ŷ
      0    1
     ---------
  0 | TN | FP |
y   | -- | -- |
  1 | FN | TP |
     ---------
"""
struct ConfusionMatrix
    TP::Int
    FN::Int
    FP::Int
    TN::Int
end

function ConfusionMatrix(ŷ::AbstractArray{Bool}, y::AbstractArray{Bool})
    ConfusionMatrix(tp(ŷ,y), fn(ŷ,y), fp(ŷ,y), tn(ŷ,y))
end
function ConfusionMatrix(ŷ::AbstractArray{<:Integer}, y::AbstractArray)
    ConfusionMatrix(Bool.(ŷ), Bool.(y))
end
function ConfusionMatrix(ŷ::AbstractArray{T}, y::AbstractArray) where T<:AbstractFloat
    ConfusionMatrix(ŷ .> one(T)/2, Bool.(y))
end

function Base.iterate(c::ConfusionMatrix, state=(c.TP, 1))
    element, current = state
    current == 2 && (element = c.FN)
    current == 3 && (element = c.FP)
    current == 4 && (element = c.TN)
    current >  4 && return nothing
    element, (element, current+1)
end

# Metrics for `ConfusionMatrix`
accuracy(c::ConfusionMatrix) = accuracy(c...)
precision(c::ConfusionMatrix) = precision(c.TP, c.FP)
recall(c::ConfusionMatrix) = recall(c.TP, c.FN)
f1(c::ConfusionMatrix) = f1(c.TP, c.FN, c.FP)
mcc(c::ConfusionMatrix) = mcc(c...)
FisherExactTest(c::ConfusionMatrix) = FisherExactTest(c...)

# Metrics for `ŷ` and `y`
for f = (:accuracy, :precision, :recall, :f1, :FisherExactTest, :mcc)
    @eval ($f)(ŷ::AbstractArray, y::AbstractArray) = ($f)(ConfusionMatrix(ŷ, y))
end

# performance summary

struct Performance
    accuracy::Float64
    precision::Float64
    recall::Float64
    f1::Float64
    mcc::Float64
    # fisher::Float64
end

function Performance(ŷ::AbstractArray{Bool}, y::AbstractArray{Bool})
    c = ConfusionMatrix(ŷ[:], y[:])
    Performance(
        accuracy(c),
        precision(c),
        recall(c),
        f1(c),
        mcc(c),
        # pvalue(FisherExactTest(c)),
        )
end
function Performance(s::AbstractArray{T}, y::AbstractArray{<:Real}) where T<:Real
    Performance(s .> one(T)/2, Bool.(y))
end
function Performance()
    Performance(0., 0., 0., 0., 0.)
end

function Base.show(io::IO, p::Performance)
    println(io, "Performance(ŷ, y)")

    for i = fieldnames(typeof(p))
        @printf(io, "%9s %.3f\n", i, getfield(p, i))
    end
end

for f = ("accuracy", "precision", "recall", "f1", "mcc", "fisher")
    @eval $(Symbol(f))(p::Performance) = getfield(p, Symbol($f))
end

# mean of a vector of `Performance`s as an instance of `Performance`
function Statistics.mean(ps::Vector{Performance})
    Performance([mean(eval(f),ps) for f = fieldnames(Performance)]...)
end

function DataFrames.DataFrame(ps::Vector{Performance})
    function f(p::Performance)
        metrics = [fieldnames(typeof(p))...]
        values = map(x->getfield(p, x), metrics)
        metrics, values
    end

    metrics_values = f.(ps)
    metrics=string.(vcat(getindex.(metrics_values, 1)...))
    values=vcat(getindex.(metrics_values, 2)...)
    DataFrame(metrics=metrics, values=values)
end

# precision-recall curves

function precision_recall(s::AbstractVector{Float64},
                          y::AbstractVector{T},
                          pos::Set{T}) where T
    n = length(s) + 1
    ordering = sortperm(s, rev=true)
    y = y[ordering]
    npos = count(x->x ∈ pos, y)
    nneg = length(y) - npos
    tn, fn, tp, fp = 0, 0, npos, nneg

    p = Array{Float64}(undef, n)
    r = Array{Float64}(undef, n)
    p[n] = precision(tp, fp)
    r[n] = recall(tp, fn)
    auc, prev_r = 0., r[n]
    pmin, pmax = p[n], p[n]

    for i in n-1:-1:1  # lowest to highest score
        dtn = y[i] ∈ pos ? 0 : 1
        dfn = 1 - dtn
        tn += dtn
        fn += dfn
        tp = npos - fn
        fp = nneg - tn

        p[i] = (tp+fp) == 0 ? dfn : precision(tp, fp)
        r[i] = recall(tp, fn)

        # Update max precision observed for current recall value
        if r[i] == prev_r
            pmax = p[i]
        else
            pmin = p[i]
            auc += (pmin + pmax) / 2 * (prev_r - r[i]) # trapezoidal rule
            prev_r = r[i]
            pmax = p[i]
        end
    end
    auc, p, r
end

"""
    PR(s, y, pos)

Calculate the precision-recall curve and its AUC using `y` labels and `s` scores for the
positive class(es) `pos`.
"""
struct PR
    auc::Float64
    precision::Vector{Float64}
    recall::Vector{Float64}
    n::Int
    p::Int
end
function PR(s::AbstractArray{Float64}, y::AbstractArray{T}, pos::Set{T}) where T
    PR(precision_recall(s[:], y[:], pos)..., length(y), sum(y))
end
function PR(s::AbstractArray{Float64}, y::AbstractArray{T}, pos::T=maximum(y)) where T
    PR(precision_recall(s[:], y[:], Set([pos]))..., length(y), sum(y))
end

for f = ("precision", "recall", "auc")
    @eval $(Symbol(f))(p::PR) = getfield(p, Symbol($f))
end

# Requires Plots.jl
macro plotpr(PRs)
    esc(quote
        local _PRs = $PRs
        if _PRs isa PR
            _PRs = [_PRs]
        end
        local _pr = _PRs[1]
        local _baseline = _pr.p / _pr.n

        fig = plot($recall(_pr), $precision(_pr), label="AUC = $(round(auc(_pr), digits=3))",
            xlabel="Recall", ylabel="Precision", size=(400,350), legend=:topright,
            xlim=(0,1), ylim=(0,1))

        for i = 2:length(_PRs)
            _pr = _PRs[i]
            plot!($recall(_pr), $precision(_pr), label="AUC = $(round(auc(_pr), digits=3))")
        end

        hline!([_baseline], label="Baseline $(round(_baseline, digits=3))", c=:black,
            ls=:dash)

        fig
    end)
end

# macro plotpr(PR)
#     esc(quote
#         plot($recall($PR), $precision($PR), label=string(round(auc($PR), digits=3)),
#             xlabel="Recall", ylabel="Precision", size=(400,350), legend=:topright)
#         baseline = $PR.p / $PR.n
#         hline!([baseline], label="Baseline $(round(baseline, digits=3))", c=:black, ls=:dash)
#     end)
# end

# Load ML

struct MLFileCollection <: AbstractFileCollection end

const ML = MLFileCollection()

function load(::MLFileCollection, T::Type=DataFrame;
              networkembeddings::Bool=false,
              center::Bool=false)
    Xs = []
    push!(Xs, load(GrowthPhenotypesWideform))
    networkembeddings && push!(Xs, load(NetworkEmbeddings))
    length(Xs) > 1 ? (X = join(Xs..., on=:id)) : (X = Xs[1])
    center && center!(X)

    Y = load(GeneOntology.GOSlimTargets)

    commonids = sort(X[:id] ∩ Y[:id])
    X = sort!(X[map(id->id ∈ commonids, X[:id]), :], :id)
    Y = sort!(Y[map(id->id ∈ commonids, Y[:id]), :], :id)

    if T <: AbstractMatrix
        return permutedims.(convert.(T{Float64}, (X[2:end], Y[2:end])))..., names(Y)[2:end]
    end

    X, Y
end

# Gene network embeddings
@file NetworkEmbeddings "$(ENV["POMBAGEDB"])/Data/network_embeddings/network_embeddings.csv"

load(f::NetworkEmbeddingsFile) = DataFrame(load(filepath(f)))
