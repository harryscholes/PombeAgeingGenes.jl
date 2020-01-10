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
            # score = parse(Float64, l[13]) # independent E-value
            score = parse(Float64, l[14]) # bit score
            # if score < threshold
            #     push!(hits, FunFamHit(l[1][1:end-4], l[4], score))
            # end
            push!(hits, FunFamHit(l[1][1:end-4], l[4], score))
        end
    end
    return hits
end

"Load FunFam hits from a set of FunFams."
function load(::FunFamHitsFile, ffs::Union{AbstractVector,AbstractSet};
              threshold::AbstractFloat=0.001)
    ffs = Set(ffs)
    hits = load(FunFamHits; threshold=threshold)
    return hits[map(x->x.ff in ffs, hits)]
end

"Load FunFam hits for all FunFams that have proteins annotated with a GO term."
function load(x::FunFamHitsFile, goterm::AbstractString; threshold::AbstractFloat=0.001)
    ffs = vec(readdlm(joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
        "funfams_with_goterms", goterm*".csv"), String))
    return load(FunFamHits, ffs; threshold=threshold)
end

# Convert `AbstractVector{FunFamHit}` into a wideform DataFrame of ids × FunFams
function Base.convert(::Type{T}, hits::AbstractVector{FunFamHit}) where T<:AbstractDataFrame
    df = T(hits)
    unique!(df, [:id, :ff])
    df = unstack(df, :id, :ff, :score)
    for col = names(df)
        col == :id && continue
        df[!, col] = -log10.(df[:, col])
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
