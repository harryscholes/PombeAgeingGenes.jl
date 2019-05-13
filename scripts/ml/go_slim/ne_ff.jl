#=
Predict GO Slim terms using:
    - network embeddings
    - FunFams
=#

include("src.jl")

dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/go_slim/RandomForestClassifier/ne_ff"
isdir(dir) || mkpath(dir)

X, Y, goterms = load(ML, Matrix, growthphenotypes=false, networkembeddings=true)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

function cvgoslim(X, Y, goterms)
    for i = 1:size(Y,1)
        goterm = goterms[i]
        try
            X, Y, goterms = load(ML, Matrix, growthphenotypes=false, networkembeddings=true,
                funfam=string(goterm))
        # An ArgumentError is thrown when the funfams_with_goterms file is empty and no
        # FunFams have proteins with a particular GO term.
        catch ArgumentError
            println("No FunFams have GO term $goterm")
            continue
        end
        cvgoslim(X, Y, goterms, i)
    end
end

cvgoslim(X, Y, goterms)
