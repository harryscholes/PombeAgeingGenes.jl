#=
Predict GO Slim terms using:
    - network embeddings
    - FunFams

session32
cd git/PombeAgeingGenes.jl/scripts/ml/go_slim/prediction
export JULIA_NUM_THREADS=32
export N_TREES=100
julia --project=$GIT/pombage
=#

include("src.jl")

dir = setupdir("prediction/gp_ne_cor")
cvpredictgoterms(dir=dir, cor=true)

dir = setupdir("prediction/gp_ne_cor2")
cvpredictgoterms(dir=dir, cor2=true)
