module PombeAgeingGenes

using DataFrames, Statistics, StatsBase, MLBase, Printf, HypothesisTests,
    Random, CSVFiles, MLDataUtils, DecisionTree, Distances, Clustering, Plots,
    Plots.PlotMeasures # TODO FIX ME

import FileIO: load
import Base: precision
import MLDataUtils: default_obsdim
import Clustering: nnodes
import StatsPlots: treepositions

include("utils.jl")

export
    @in

include("io.jl")

export
    AbstractFile,
    AbstractFileCollection,
    @file,
    filepath,
    load

include("geneontology.jl")

export
    GeneOntology

include("stats.jl")

export
    snr,
    coefvar

include("ml.jl")

export
    traintestsplit,
    standardize,
    softmax,
    CrossValidation,
    fit,
    ML,
    stratifiedkfolds,
    tp,
    fn,
    fp,
    tn,
    accuracy,
    precision,
    recall,
    f1,
    Performance,
    PR,
    @plotpr,
    mcc,
    auc

include("fisherlineardiscriminant.jl")

export
    FisherLinearDiscriminant,
    evaluate,
    probability,
    predict

include("growthphenotypes.jl")

export
    controls,
    GrowthPhenotypes,
    GrowthPhenotypesNoOutliers,
    GrowthPhenotypesWideform,
    meansizes,
    impute!,
    isoutlier,
    findoutliers!,
    findoutliers,
    outliers!,
    outliers,
    removeoutliers!,
    removeoutliers,
    findrepeats!,
    findrepeats,
    removerepeats!,
    removerepeats

include("plotting.jl")

export
    clusteredheatmap

end # module
