include("../src.jl")

macro crossvalidatepredict(args...)
    X = args[1]
    y = args[2]
    ids = args[3]
    ex = args[end]

    k = 5
    f = stratifiedkfolds

    esc(quote
        local data = shuffleobs(($X, $ids, $y))
        local ys = []
        local ŷs = []
        local yids = []

        for ((Xtrain, _, ytrain), (Xtest, idtest, ytest)) = ($f)(data, $k)
            Xtrain, ytrain, Xtest, ytest, idtest =
                map(collect, (Xtrain, ytrain, Xtest, ytest, idtest))
            ŷ = $ex
            push!(ys, ytest)
            push!(ŷs, ŷ)
            push!(yids, idtest)
        end

        ŷs, ys, yids, Performance.(ŷs, ys), PR(vcat(ŷs...), vcat(ys...))
    end)
end

function cvpredictgoterms(; dir, runnumber=0)
    goterms = map(x->replace(x, ":"=>""), GeneOntology.GO_SLIM_TERMS)

    for (i, goterm) = enumerate(goterms)
        println("goterm = ", goterm)

        X, Y = try
            load(ML, growthphenotypes=false, networkembeddings=true, funfam=goterm)
        catch ArgumentError
            @warn "Error"
            continue
        end

        ids = X[:id]
        X = permutedims(Matrix(X[2:end]))
        y = [j == 1 for j = Y[Symbol(goterm)]]

        grid = Dict(
            :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
            :partial_sampling => [.5, .75, 1.],
            )

        Random.seed!(runnumber+i)
        cvpredictgoterm(X, y, grid, ids, goterm; dir=dir)
    end
end

function cvpredictgoterm(X, y, grid, ids, goterm; dir)
    fname = joinpath(dir, goterm)
    ŷs, ys, yids, p, pr = cvpredict(X, y, ids, grid)
    savefig(plot(pr, legend=:topright), fname*".pdf")
    writecvresults(fname*".json", ŷs, ys, "yids"=>yids)
end

function cvpredict(X, y, ids, grid)
    ŷs, ys, yids, p, pr = @crossvalidatepredict X y ids begin
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        params = @gridsearch rfc grid
        ŷ = rfc(Xtrain, ytrain, Xtest; params...)
    end
end
