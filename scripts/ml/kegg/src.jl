using PombeAgeingGenes, Distributed, DecisionTree, JSON, Plots

using Distributed, DecisionTree

if haskey(ENV, "NPROCS")
    addprocs(parse(Int, ENV["NPROCS"]))
    @show nprocs()
end

@everywhere using DecisionTree

commit_hash = "8e62edf"

basedir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "kegg",
                   "RandomForestClassifier", commit_hash)

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

function cvkegg(X, Y, pathways)
    for i = 1:size(Y,1)
        cvkegg(X, Y, pathways, i)
    end
end

function cvkegg(X, Y, pathways, i)
    pathway = pathways[i]
    @show pathway
    fname = "$(dir)/$pathway"
    y = [i == 1 for i = Y[i,:]]
    count(y) < 5 && return nothing
    ŷs, ys, p, pr = cv(X, y, grid)
    savefig(plot(pr, legend=:topright), fname * ".pdf")
    writecvresults(fname * ".json", ŷs, ys)
end
