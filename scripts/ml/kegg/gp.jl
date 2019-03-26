#=
Predict KEGG pathways using:
    - growth phenotypes
=#

include("src.jl")

dir = joinpath(basedir, "gp")
isdir(dir) || mkpath(dir)

X, Y, pathways = load(ML, Matrix, growthphenotypes=true, networkembeddings=false, Y=:kegg)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

cvkegg(X, Y, pathways)
