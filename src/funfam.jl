struct FunFamHit
    id::String
    ff::String
    score::Float64
end

@file FunFamHits joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam", "peptide.domtblout")

"Load all FunFam hits."
function load(T::FunFamHitsFile; threshold::AbstractFloat=0.001)
    hits = FunFamHit[]

    open(filepath(T), "r") do io
        for line = eachline(io)
            (startswith(line, "#") || length(line) == 0) && continue
            l = split(line)
            score = parse(Float64, l[13])

            if score < threshold
                push!(hits, FunFamHit(l[1][1:end-4], l[4], score))
            end
        end
    end

    hits
end

"Load FunFam hits from a set of FunFams."
function load(::FunFamHitsFile, ffs::Union{AbstractVector,AbstractSet};
              threshold::AbstractFloat=0.001)
    ffs = Set(ffs)
    hits = load(FunFamHits; threshold=threshold)
    hits[map(x->x.ff in ffs, hits)]
end

"Load FunFam hits for all FunFams that have proteins annotated with a GO term."
function load(x::FunFamHitsFile, goterm::AbstractString; threshold::AbstractFloat=0.001)
    ffs = vec(readdlm(joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
        "funfams_with_goterms", goterm*".csv"), String))
    load(FunFamHits, ffs; threshold=threshold)
end

# Convert `AbstractVector{FunFamHit}` into a wideform DataFrame of ids × FunFams
function Base.convert(::Type{T}, hits::AbstractVector{FunFamHit}) where T<:AbstractDataFrame
    df = T(hits)
    unique!(df, [:id, :ff])
    df = unstack(df, :id, :ff, :score)

    for col = names(df)
        col == :id && continue
        df[col] = -log10.(df[col])
        df[col] = coalesce.(df[col], 0.)
    end

    df
end
