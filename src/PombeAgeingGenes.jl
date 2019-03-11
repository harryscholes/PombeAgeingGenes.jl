module PombeAgeingGenes

using
    DataFrames,
    Statistics,
    StatsBase,
    MLBase,
    Printf,
    Random,
    CSVFiles,
    MLDataUtils

import FileIO: load
import Base: precision
import MLDataUtils: default_obsdim

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
    stratifiedkfolds,
    @crossvalidate,
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
    auc,
    ML,
    NetworkEmbeddings

# include("fisherlineardiscriminant.jl")
#
# export
#     FisherLinearDiscriminant,
#     evaluate,
#     probability,
#     predict

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

# using Distances, Clustering, Plots, Plots.PlotMeasures
#
# import Clustering: nnodes
# import StatsPlots: treepositions
#
# include("plotting.jl")
#
# export
#     clusteredheatmap

end # module
