#=
Predict GO Slim terms using:
    - growth phenotypes
    - network embeddings
=#

include("src.jl")

dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/go_slim/RandomForestClassifier/gp_ne"
isdir(dir) || mkpath(dir)

X, Y, goterms = load(ML, Matrix, networkembeddings=true)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

cvgoslim(X, Y, goterms)
