#=
Predict GO Slim terms using:
    - network embeddings
    - FunFams
=#

include("src.jl")

dir = setupdir("ne_ff")

X, Y, goterms = load(ML, Matrix, growthphenotypes=false, networkembeddings=true)

const grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .75, 1.],
    )

function cvgoterms(X, Y, goterms; dir, runnumber=0)
    for i = 1:size(Y,1)
        goterm = string(goterms[i])
        @show goterm
        try
            X, Y, goterms = load(ML, Matrix, growthphenotypes=false,
                                 networkembeddings=true, funfam=goterm)
        # An ArgumentError is thrown when the funfams_with_goterms file is empty and no
        # FunFams have proteins with a particular GO term.
        catch ArgumentError
            println("No FunFams have GO term $goterm")
            continue
        end
        y = [j == 1 for j = Y[i,:]]
        Random.seed!(runnumber+i)
        cvgoterm(X, y, goterm; dir=dir)
    end
end

cvgoterms(X, Y, goterms; dir=dir)
repeats(X, Y, goterms, 5; dir=dir)
