using PombeAgeingGenes, Distributed, DecisionTree, JSON, Plots, Random

using Distributed, DecisionTree

if haskey(ENV, "NPROCS")
    addprocs(parse(Int, ENV["NPROCS"]))
    @show nprocs()
end

@everywhere using DecisionTree

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

function cvgoslim(X, Y, goterms)
    for i = 1:size(Y,1)
        cvgoslim(X, Y, goterms, i)
    end
end

function cvgoslim(X, Y, goterms, i)
    goterm = goterms[i]
    @show goterm
    fname = "$dir/$(goterm)"
    y = [i == 1 for i = Y[i,:]]
    Random.seed!(i)
    ŷs, ys, p, pr = cv(X, y, grid)
    savefig(plot(pr, legend=:topright), fname * ".pdf")
    writecvresults(fname * ".json", ŷs, ys)
end
