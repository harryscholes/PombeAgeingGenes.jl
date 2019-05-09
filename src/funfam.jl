struct FunFamHit
    id::String
    ff::String
    score::Float64
end

@file FunFamHits joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam", "peptide.domtblout")

"""Load all FunFam hits."""
function load(x::FunFamHitsFile; threshold::AbstractFloat=0.001)
    hits = FunFamHit[]
    open(filepath(x), "r") do io
        for line = eachline(io)
            startswith(line, "#") && continue
            length(line) == 0 && continue
            l = split(line)
            score = parse(Float64, l[13])
            if score < threshold
                push!(hits, FunFamHit(l[1][1:end-4], l[4], score))
            end
        end
    end
    df = finalize(hits)
end

"""Load FunFam hits from a set of FunFams."""
function load(x::FunFamHitsFile, ffs::Union{AbstractVector,AbstractSet};
              threshold::AbstractFloat=0.001)
    ffs = Set(ffs)
    hits = FunFamHit[]
    open(filepath(x), "r") do io
        for line = eachline(io)
            startswith(line, "#") && continue
            length(line) == 0 && continue
            l = split(line)
            ff = l[4]
            if ff in ffs
                score = parse(Float64, l[13])
                if score < threshold
                    push!(hits, FunFamHit(l[1][1:end-4], ff, score))
                end
            end
        end
    end
    df = finalize(hits)
end

"""Load FunFam hits for all FunFams that have proteins annotated with a GO term."""
function load(x::FunFamHitsFile, goterm::AbstractString; threshold::AbstractFloat=0.001)
    ffs = vec(readdlm(joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
        "funfams_with_goterms", goterm*".csv"), String))
    load(x, ffs; threshold=threshold)
end

function Base.finalize(hits::AbstractVector{FunFamHit})
    df = DataFrame(hits)
    unique!(df, [:id, :ff])
    df = unstack(df, :id, :ff, :score)
    for col = names(df)
        col == :id && continue
        df[col] = -log10.(df[col])
        df[col] = coalesce.(df[col], 0.)
    end
    df
end
