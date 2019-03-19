#=
Predict GO Slim terms using:
    - growth phenotypes
=#

include("src.jl")

dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/go_slim/RandomForestClassifier/gp"
isdir(dir) || mkpath(dir)

X, Y, goterms = load(ML, Matrix)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

cvgoslim(X, Y, goterms)