#=
Predict GO Slim terms using:
    - network embeddings
    - FunFams
=#

include("src.jl")

dir = setupdir("prediction/ne_ff")

cvpredictgoterms(dir=dir)
