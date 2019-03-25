#=
S. pombe KEGG pathways
=#

@file KEGGPathways joinpath(ENV["POMBEAGEINGGENES"], "data", "kegg", "kegg_pathways.tsv")

function load(x::KEGGPathwaysFile)
    d = Dict{String,String}()
    open(filepath(x), "r") do io
        for line = eachline(io)
            line = split(line, '\t')
            d[line[1][6:end]] = line[2]
        end
    end
    d
end

const KEGGPATHWAYS = ["spo00010", "spo00020", "spo00030", "spo00040", "spo00051",
    "spo00052", "spo00061", "spo00062", "spo00071", "spo00072", "spo00100", "spo00130",
    "spo00190", "spo00220", "spo00230", "spo00240", "spo00250", "spo00260", "spo00261",
    "spo00270", "spo00280", "spo00290", "spo00300", "spo00310", "spo00330", "spo00332",
    "spo00340", "spo00350", "spo00360", "spo00380", "spo00400", "spo00410", "spo00430",
    "spo00440", "spo00450", "spo00460", "spo00472", "spo00480", "spo00500", "spo00510",
    "spo00511", "spo00513", "spo00514", "spo00515", "spo00520", "spo00561", "spo00562",
    "spo00563", "spo00564", "spo00565", "spo00590", "spo00600", "spo00620", "spo00630",
    "spo00640", "spo00650", "spo00660", "spo00670", "spo00680", "spo00730", "spo00740",
    "spo00750", "spo00760", "spo00770", "spo00780", "spo00785", "spo00790", "spo00860",
    "spo00900", "spo00909", "spo00910", "spo00920", "spo00970", "spo01040", "spo01100",
    "spo01110", "spo01130", "spo01200", "spo01210", "spo01212", "spo01230", "spo01502",
    "spo02010", "spo03008", "spo03010", "spo03013", "spo03015", "spo03018", "spo03020",
    "spo03022", "spo03030", "spo03040", "spo03050", "spo03060", "spo03410", "spo03420",
    "spo03430", "spo03440", "spo03450", "spo04011", "spo04070", "spo04111", "spo04113",
    "spo04120", "spo04122", "spo04130", "spo04136", "spo04138", "spo04139", "spo04141",
    "spo04144", "spo04145", "spo04146", "spo04933"]

struct KEGGPathway <: AbstractFile
    pathway::String
end

load(x::KEGGPathway) =
    JSON.parsefile(joinpath(ENV["POMBEAGEINGGENES"], "data", "kegg", x.pathway*".json"))[1]

"""
    kegg"<pathway id>"
"""
macro kegg_str(pathway)
    :(KEGGPathway($pathway))
end

@file KEGGPathwayTargets joinpath(ENV["POMBEAGEINGGENES"], "data",
                                   "keggpathway_targets.csv")

load(x::KEGGPathwayTargetsFile) = DataFrame(load(filepath(x)))
