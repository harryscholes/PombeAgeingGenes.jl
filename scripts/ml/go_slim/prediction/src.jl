# Without gridsearch

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

function cvpredictgoterms(; dir, runnumber=0, kwargs...)
    goterms = string.(load(ML, Matrix;
        growthphenotypes=true,
        networkembeddings=true,
        kwargs...
        )[3])

    for (i, goterm) = enumerate(goterms)
        isfile(joinpath(dir, goterm*".json")) && continue # skip predictions that have already been made

        println("goterm = ", goterm)

        X, Y = load(ML;
            growthphenotypes=true,
            networkembeddings=true,
            funfam=replace(goterm, "GO"=>"GO:"),
            kwargs...
            )

        ids = X[:,:id]
        X = permutedims(Matrix(X[:,2:end]))
        y = [j == 1 for j = Y[:,Symbol(goterm)]]

        all(.!y) && continue # no genes annotated with go term. Why does this happen?

        Random.seed!(runnumber+i)
        cvpredictgoterm(X, y, ids, goterm; dir=dir)
    end
end

function cvpredictgoterm(X, y, ids, goterm; dir)
    fname = joinpath(dir, goterm)
    ŷs, ys, yids, p, pr = cvpredict(X, y, ids)
    savefig(plot(pr, legend=:topright), fname*".pdf")
    writecvresults(fname*".json", ŷs, ys, "yids"=>yids)
end

function cvpredict(X, y, ids)
    ŷs, ys, yids, p, pr = @crossvalidatepredict X y ids begin
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        ŷ = rfc(Xtrain, ytrain, Xtest)
    end
end

#=
# With gridsearch

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
    goterms = string.(load(ML, Matrix, growthphenotypes=false, networkembeddings=true)[3])
    for (i, goterm) = enumerate(goterms)
        isfile(joinpath(dir, goterm*".json")) && continue # skip predictions that have already been made

        println("goterm = ", goterm)

        X, Y = load(ML, growthphenotypes=false, networkembeddings=true, funfam=replace(goterm, "GO"=>"GO:"))
        ids = X[:,:id]
        X = permutedims(Matrix(X[:,2:end]))
        y = [j == 1 for j = Y[:,Symbol(goterm)]]
        all(.!y) && continue # no genes annotated with go term. Why does this happen?
        X_sqrt = floor(sqrt(size(X,1)))
        n_subfeatures = [10, 25, 50]
        X_sqrt < 50 && push!(n_subfeatures, X_sqrt)
        grid = Dict(
            :n_subfeatures => n_subfeatures,
            :partial_sampling => [.5, .75, 1.],
            )
        ŷ = rfc(Xtrain, ytrain, Xtest; params...)
    end
end
=#
