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
    Requires,
    DelimitedFiles,
    CodecZlib,
    DataStructures,
    OBOParse

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
    PValue,
    Greater,
    Less

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
    # GrowthPhenotypesTrigitised,
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
    removerepeats,
    trigitise!

include("kegg.jl")

export
    KEGGPathways,
    KEGGPATHWAYS,
    KEGGPathway,
    @kegg_str,
    KEGGPathwayTargets

include("funfam.jl")

export
    FunFamHits,
    FunFamGOTerms,
    FunFamGOTermsIEA,
    FunFamInclusionThresholds

end # module
