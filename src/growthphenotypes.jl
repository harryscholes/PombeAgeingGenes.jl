const controls = ["wt", "SPBC29B5.01", "SPBC106.10", "SPBC1105.14", "SPAC1687.15",
    "SPBC725.11c", "SPCC18B5.11c", "SPBC1198.14c"]

# IO

#=
Jan2019 growth phenotypes
=#

#=
const _fname = joinpath(ENV["POMBEAGEINGGENES"], "data", "Jan2019_BBSRC_results")

# Created by $POMBEAGEINGGENES/scripts/growth_phenotypes/process.jl
@file(GrowthPhenotypes, _fname * "_clean.csv")
@file(GrowthPhenotypesNoOutliers, _fname * "_no_outliers.csv")
@file(GrowthPhenotypesWideform, _fname * "_no_outliers_wideform.csv")
@file(GrowthPhenotypesTrigitised, _fname * "_no_outliers_wideform_trigitised.csv")

for T = (GrowthPhenotypesFile, GrowthPhenotypesNoOutliersFile)
    @eval function load(x::$T)
        df = DataFrame(load(filepath(x), colparsers=Dict(:hash=>UInt64)))
        df[!, :phlox] = parse.(Bool, df[:, :phlox])
        return df
    end
end

# load(x::GrowthPhenotypesWideformFile) = DataFrame(load(filepath(x)))

for T = (GrowthPhenotypesWideformFile, GrowthPhenotypesTrigitisedFile)
    @eval load(x::$T) = DataFrame(load(filepath(x)))
end
=#

#=
Oct2019 growth phenotypes
=#

const _fname = joinpath(ENV["POMBEAGEINGGENES"], "data", "Oct2019_BBSRC_results")

# Created by $POMBEAGEINGGENES/scripts/growth_phenotypes/process_2.jl
@file(GrowthPhenotypes, _fname * "_clean.csv")
@file(GrowthPhenotypesNoOutliers, _fname * "_normalised.csv")
@file(GrowthPhenotypesWideform, _fname * "_wideform.csv")

function load(x::GrowthPhenotypesFile)
    df = DataFrame(load(filepath(x)))
    df.phlox = parse.(Bool, df.phlox)
    return df
end

for T = (GrowthPhenotypesNoOutliersFile, GrowthPhenotypesWideformFile)
    @eval load(x::$T) = DataFrame(load(filepath(x)))
end

# Growth phenotypes processing

"""
    meansizes(df[; nrepeats, digits])

Calculate mean sizes per strain-condition pair.

If the number of repeats is ≤ `nrepeats` then the size is set to `missing`. Sizes are rounded to
`digits` digits.
"""
function meansizes(df::DataFrame; nrepeats::Int=2, digits::Int=2)
    return by(df, [:id, :condition],
        size = :size => x->length(x) ≤ nrepeats ? missing : round(mean(x), digits=digits))
end

"""
    impute!(df)

Impute `NaN`s with mean size per condition.
"""
function impute!(df::DataFrame)
    for g = groupby(df, :condition)
        g[ismissing.(g[!,:size]), :size] .= mean(skipmissing(g[:, :size]))
    end
    return df
end

# QC

"""
    Fence(low, high)
    Fence(xs[; scale])

Outlier fence.
"""
struct Fence{T<:Real}
    low::T
    high::T
end
function Fence(xs::AbstractVector; scale::Real=3)
    m = median(xs)
    mad = StatsBase.mad(xs, center=m, normalize=false)
    threshold = mad * scale
    return Fence(m - threshold, m + threshold)
end

isoutlier(F::Fence{T}, x::T) where T<:Real = !(F.low ≤ x ≤ F.high)

"""
    findoutliers!(df[; scale])
    findoutliers(df[; scale])

Add a column `:isoutlier` to the DataFrame `df` that identitfies whether a colony size is an outlier for in a strain-condition pair.

Outliers are defined as being MAD * `scale` larger/smaller than then median.
"""
function findoutliers!(df::DataFrame; scale::Real=3)
    df[!, :isoutlier] .= false
    for g = groupby(df, [:id, :condition])
        f = Fence(g[!, :size]; scale=scale)
        g[:, :isoutlier] = map(x->isoutlier(f, x), g[:, :size])
    end
    return df
end
findoutliers(df::DataFrame; scale::Real=3) = findoutliers!(deepcopy(df); scale=scale)

"""
    outliers!(df; keepoutliers[, scale])
    outliers(df; keepoutliers[, scale])

Return the DataFrame `df` containing only the outliers if `keepoutliers` is `true` and
removes outliers if `keepoutliers` is false.

`scale` is passed to `findoutliers!`.
"""
function outliers!(df::DataFrame; keepoutliers::Bool, scale::Real=3)
    findoutliers!(df)
    deleterows!(df, df[!, :isoutlier] .!= keepoutliers)
    select!(df, Not(:isoutlier))
end
function outliers(df::DataFrame; keepoutliers::Bool, scale::Real=3)
    outliers!(deepcopy(df); keepoutliers=keepoutliers, scale=scale)
end

"""
    removeoutliers!(df[; scale])
    removeoutliers(df[; scale])

Return the DataFrame `df` with outliers removed.

`scale` is passed to `findoutliers!`.
"""
function removeoutliers!(df::DataFrame; scale::Real=3)
    outliers!(df; keepoutliers=false, scale=scale)
end
function removeoutliers(df; scale::Real=3)
    removeoutliers!(deepcopy(df); scale=scale)
end

"""
    findrepeats!(df)

Add a column `:nrepeats` to the DataFrame `df` containing the number of repeats for a
strain-condition pair.
"""
function findrepeats!(df::DataFrame)
    df[!, :nrepeats] .= 0
    for g = groupby(df, [:id, :condition])
        g[:, :nrepeats] .= size(g, 1)
    end
    return df
end

findrepeats(df::DataFrame) = findrepeats!(deepcopy(df))

"""
    removerepeats!(df; nrepeats)
    removerepeats(df; nrepeats)

Remove condition-strain pairs from DataFrame `df` with `nrepeats` or fewer repeats.
"""
function removerepeats!(df::DataFrame; nrepeats::Int)
    findrepeats!(df)
    deleterows!(df, df[:, :nrepeats] .≤ nrepeats)
    deletecols!(df, :nrepeats)
end
function removerepeats(df::DataFrame; nrepeats::Int)
    removerepeats!(deepcopy(df); nrepeats=nrepeats)
end

"""
    trigitise!(A; ratio)

Encode log2 growth pheontype data as +1, 0 and -1 for resistant, no phenotype and sensitive
strains, respectively.
"""
function trigitise!(A::AbstractArray; ratio::AbstractFloat=0.2)
    r = 1-ratio
    l = log2(r)
    u = log2(inv(r))
    A[A .< l] .= -1
    A[A .> u] .= 1
    A[l .≤ A .≤ u] .= 0
    return A
end

function trigitise!(df::AbstractDataFrame; ratio::AbstractFloat=0.2)
    for col = names(df)
        if eltype(df[:, col]) <: Real
            df[!, col] = trigitise!(df[:, col]; ratio=ratio)
        end
    end
    return df
end
