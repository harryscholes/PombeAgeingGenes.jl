#=
Predict GO Slim terms using:
    - growth phenotypes
    - network embeddings
=#

include("src.jl")

dir = setupdir("gp_ne")

X, Y, goterms = load(ML, Matrix, growthphenotypes=true, networkembeddings=true)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

cvgoterms(X, Y, goterms; dir=dir)
repeats(X, Y, goterms, 5; dir=dir)
