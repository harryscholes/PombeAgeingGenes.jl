using PombeAgeingGenes, Distributed, DecisionTree, JSON, Plots, Random

using Distributed, DecisionTree

if haskey(ENV, "NPROCS")
    addprocs(parse(Int, ENV["NPROCS"]))
    @show nprocs()
end

@everywhere using DecisionTree

function setupdir(features)
    dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "go_slim",
                   "RandomForestClassifier", features)
    isdir(dir) || mkpath(dir)
    dir
end

function rfc(Xtrain, ytrain, Xtest; kwargs...)
    model = DecisionTree.fit!(RandomForestClassifier(; n_trees=100, kwargs...), Xtrain, ytrain)
    ŷ = DecisionTree.predict_proba(model, Xtest)[:,2]
end

function cv(X, y, grid)
    ŷs, ys, p, pr = @crossvalidate X y begin
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        params = @gridsearch rfc grid
        ŷ = rfc(Xtrain, ytrain, Xtest; params...)
    end
end

function cvgoterms(X, Y, goterms; dir, runnumber=0)
    for i = 1:size(Y,1)
        goterm = string(goterms[i])
        @show goterm
        y = [j == 1 for j = Y[i,:]]
        Random.seed!(runnumber+i)
        cvgoterm(X, y, goterm; dir=dir)
    end
end

function cvgoterm(X, y, goterm; dir)
    fname = joinpath(dir, goterm)
    ŷs, ys, p, pr = cv(X, y, grid)
    savefig(plot(pr, legend=:topright), fname*".pdf")
    writecvresults(fname*".json", ŷs, ys)
end

repeats(nrepeats::Integer; dir::AbstractString) = repeats(1, nrepeats; dir=dir)

function repeats(start::Integer, stop::Integer; dir::AbstractString)
    for i = start:stop
        dir_i = joinpath(dir, string(i))
        isdir(dir_i) || mkpath(dir_i)
        cvgoterms(X, Y, goterms, dir=dir_i; runnumber=i)
    end
end
