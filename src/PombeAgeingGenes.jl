module PombeAgeingGenes

using
    DataFrames,
    Statistics,
    StatsBase,
    MLBase,
    Printf,
    Random,
    CSVFiles,
    MLDataUtils,
    RecipesBase,
    JSON,
    Requires

import FileIO: load
import Base: precision
import MLDataUtils: default_obsdim

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
end

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
    coefvar,
    NullDistribution,
    PValue

include("ml.jl")

export
    stratifiedkfolds,
    @crossvalidate,
    writecvresults,
    loadcvresults,
    @gridsearch,
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

end # module
