struct FunFamHit
    id::String
    ff::String
    score::Float64
end

# @file FunFamHits joinpath(ENV["POMBEAGEINGGENES"], "data", "cafa4", "peptide.domtbl.gz")
@file FunFamHits joinpath(ENV["POMBEAGEINGGENES"], "data", "cafa4", "peptide.cut_tc.domtbl.gz")


"Load FunFam hits for which `f` evaluates to `true` and `score` < `threshold`."
function load(f::Function, T::FunFamHitsFile; threshold::AbstractFloat=0.0001)
    hits = FunFamHit[]
    open(GzipDecompressorStream, filepath(T)) do io
        for line = eachline(io)
            (startswith(line, "#") || length(line) == 0) && continue
            l = split(line)

            # `f` is a function that is applied to each line in the `FunFamHitsFile`
            if f(l)
                id = replace(l[1], ":pep"=>"")
                ff = l[4]

                # independent E-value
                score = parse(Float64, l[13])
                if score < threshold
                    push!(hits, FunFamHit(id, ff, score))
                end

                # bit score
                # score = parse(Float64, l[14])
                # push!(hits, FunFamHit(id, ff, score))
            end
        end
    end
    return hits
end

"Load all FunFam hits."
function load(T::FunFamHitsFile; threshold::AbstractFloat=0.0001)
    return load(l->true, T, threshold=threshold)
end

"Load FunFam hits from a set of FunFams."
function load(T::FunFamHitsFile, ffs::Union{AbstractVector,AbstractSet};
              threshold::AbstractFloat=0.0001)
    ffs = Set(ffs)
    # FunFam ID is in `l[4]`
    return load(l->l[4] in ffs, T, threshold=threshold)
end

"Load FunFam hits for all FunFams that have proteins annotated with a GO term."
function load(x::FunFamHitsFile, goterm::AbstractString; threshold::AbstractFloat=0.0001)
    for line in eachline(joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam", "go2funfam.csv"))
        if startswith(line, goterm)
            ffs = string.(split(split(line, limit=2)[2], ";"))
            return load(FunFamHits, ffs; threshold=threshold)
        end
    end
    return FunFamHit[]
end

function map_goterms_to_funfams()
    d = DefaultDict{String, Vector{String}}([])
    for line in eachline(joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam", "go2funfam.csv"))
        go, ffs = split(line, limit=2)
        d[go] = split(ffs, ";")
    end
    return d
end

# Convert `AbstractVector{FunFamHit}` into a wideform DataFrame of ids × FunFams
function Base.convert(::Type{T}, hits::AbstractVector{FunFamHit}) where T<:AbstractDataFrame
    if length(hits) == 0
        return DataFrame(id=String[])
    end
    df = T(hits)
    unique!(df, [:id, :ff])
    # deleterows!(df, df.score .< 0) # bit scores
    df = unstack(df, :id, :ff, :score)
    for col = names(df)
        col == :id && continue
        df[!, col] = -log10.(df[:, col]) # E-values
        df[!, col] = coalesce.(df[:, col], 0.)
    end
    return df
end

@file FunFamGOTerms joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
                             "funfam_uniprot_goterms_parsed.csv")

function load(x::FunFamGOTermsFile)
    associations = Tuple{String,String}[]
    for line in eachline(filepath(x))
        push!(associations, Tuple(split(line, ",")))
    end
    return associations
end

@file FunFamGOTermsIEA joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
                                "funfam_uniprot_goterms_parsed_uniprotkb_kw_iea.csv")

function load(x::FunFamGOTermsIEAFile)
    associations = Tuple{String,String}[]
    for line in eachline(filepath(x))
        push!(associations, Tuple(split(line, ",")))
    end
    return associations
end

struct InclusionThreshold
	id::String
	tc::Float64
end

@file FunFamInclusionThresholds joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam", "funfam-hmm3-v4_2_0.lib.gz")

function load(x::FunFamInclusionThresholdsFile)
	xs = InclusionThreshold[]
	id = ""
	open(GzipDecompressorStream, filepath(x)) do io
		for line in eachline(io)
			if startswith(line, r"^NAME")
				id = split(line)[2]
			elseif startswith(line, r"^TC")
				tc = parse(Float64, split(line)[2])
				push!(xs, InclusionThreshold(id, tc))
			end
		end
	end
	return xs
end
