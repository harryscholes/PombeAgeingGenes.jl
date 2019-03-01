@recipe function f(A::AbstractMatrix, hci::Hclust, hcj::Hclust)
    seriestype := :heatmap
    A[hci.order, hcj.order]
end

function _clusteredheatmapsymmetric(A; metric, linkage, cor=false, kwargs...)
    hc = hclust(pairwise(metric, A, dims=2), linkage=linkage)
    plot(A, hc, hc;
        c = cor ? :RdBu_r : nothing, clims = cor ? (-1., 1.) : nothing,
        kwargs...)
end

function _clusteredheatmapasymmetric(A; metric, linkage, kwargs...)
    hci, hcj = map(d->hclust(pairwise(metric, A, dims=d), linkage=linkage), (1,2))
    plot(A, hci, hcj; kwargs...)
end

function clusteredheatmap(f, A; metric=Euclidean(), linkage=:average, kwargs...)
    A = f(A)
    if issymmetric(A)
        _clusteredheatmapsymmetric(A; metric=metric, linkage=linkage, cor=f==cor, kwargs...)
    else
        _clusteredheatmapasymmetric(A; metric=metric, linkage=linkage, kwargs...)
    end
end

function clusteredheatmap(A; metric=Euclidean(), linkage=:average, kwargs...)
    clusteredheatmap(identity, A; metric=metric, linkage=linkage, kwargs...)
end
