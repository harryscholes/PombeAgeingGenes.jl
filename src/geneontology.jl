module GeneOntology

export
    RELATIONSHIPS,
    EVIDENCE_CODES,
    EVIDENCE_CODES_TRUSTED,
    GO,
    GOAnnotations,
    GO_SLIM_TERMS,
    slim2descendants,
    descendants2slim,
    GOSlimTargets,
    AGEING_FYPO_TERMS

using PombeAgeingGenes, OBOParse, DataFrames

import PombeAgeingGenes: load

# Ontology
const GO_RELATIONSHIPS = [:part_of, :is_a, :regulates, :positively_regulates,
   :negatively_regulates]

const EVIDENCE_CODES = Dict(
   :experimental => ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"],
   :curated      => ["ISS", "ISO", "ISA", "ISM", "IGC", "IBA",
                     "IBD", "IKR", "IRD", "RCA", "TAS", "IC"],
   :automatic    => ["IEA"])

const EVIDENCE_CODES_TRUSTED = [EVIDENCE_CODES[:experimental]; EVIDENCE_CODES[:curated]]

function OBOParse.descendants(ontology::Ontology, term::AbstractString)
    map(x->x.id, descendants(ontology, ontology[term], GO_RELATIONSHIPS))
end

@file(GO, ENV["POMBEAGEINGGENES"] * "/data/go-basic.obo")

load(F::GOFile) = OBOParse.load(filepath(F), "GO")

# Annotations
@file(GOAnnotations, ENV["POMBEAGEINGGENES"] * "/data/gene_association.pombase")

struct Annotation
    id::String
    go::String
    ec::String
end

function load(F::GOAnnotationsFile; ecs=EVIDENCE_CODES_TRUSTED)
    annotations = Annotation[]

    open(F) do f
        for l = eachline(f)
            startswith(l, "!") && continue # comment lines

            s = split(l)
            a = Annotation(s[2], s[4], s[6])

            a.ec ∉ ecs               && continue
            !startswith(a.id, "SP")  && continue
            !startswith(a.go, "GO:") && continue

            push!(annotations, a)
        end
    end

    annotations
end

# Slim
const GO_SLIM_TERMS = [
   "GO:0030036", "GO:0006915", "GO:0030437", "GO:0006914", "GO:0005975", "GO:0006766",
   "GO:0007155", "GO:0071554", "GO:0006520", "GO:0006325", "GO:0051186", "GO:0000747",
   "GO:0002181", "GO:0098754", "GO:0006310", "GO:0006281", "GO:0006260", "GO:0007163",
   "GO:0006091", "GO:0006629", "GO:0140013", "GO:0061024", "GO:0055065", "GO:0000226",
   "GO:0140053", "GO:0007005", "GO:0000281", "GO:0000070", "GO:0016071", "GO:0071941",
   "GO:0055086", "GO:0006913", "GO:0140056", "GO:0007031", "GO:0030163", "GO:0006457",
   "GO:0006486", "GO:0051604", "GO:0070647", "GO:0006605", "GO:0065003", "GO:1901990",
   "GO:0006355", "GO:0042254", "GO:0023052", "GO:0016074", "GO:0016073", "GO:0006790",
   "GO:0032200", "GO:0006351", "GO:0055085", "GO:0006399", "GO:0016192"]

"""
Map GO Slim terms to their descendent terms in the GO.
"""
slim2descendants(ontology) = Dict(t=>descendants(ontology, t) for t = GO_SLIM_TERMS)

"""
Map GO terms to all GO Slim ancestor terms.
"""
function descendants2slim(s2d)
    d = Dict{String,Set{String}}()
    for (k, vs) = s2d, v = vs
        haskey(d, v) ? push!(d[v], k) : (d[v] = Set([k]))
    end
    d
end

@file(GOSlimTargets, ENV["POMBEAGEINGGENES"] * "/data/goslim_targets.csv")

function load(x::GOSlimTargetsFile)
    df = DataFrame(load(filepath(x)))

    for col = names(df)
        df[!, col] = coalesce.(df[:, col], 0)
    end

    df
end

load(F::GOSlimTargetsFile, goterm::Symbol) = load(F)[[:id, goterm]]

"""
Ageing FYPO terms

- increased viability upon nitrogen starvation (FYPO:0004344)
- increased viability in stationary phase (FYPO:0001309)
"""
const AGEING_FYPO_TERMS = Dict(
    :FYPO0004344 => [
        "SPBC1861.02", "SPBC2G2.06c", "SPBC691.03c", "SPBC685.04c", "SPBC1D7.03",
        "SPBC428.08c", "SPAC4G9.11c", "SPAC23A1.06c", "SPAC1805.07c", "SPBC947.05c",
        "SPAC17A5.16", "SPCC1753.02c", "SPBC32F12.03c", "SPAC1687.15", "SPCC4G3.09c",
        "SPAC1834.04", "SPBC947.09", "SPBC19C7.01", "SPAC694.06c", "SPAC8F11.03",
        "SPAC1F3.09", "SPAC18B11.04", "SPCC663.06c", "SPBC543.07", "SPAC1F8.06",
        "SPBC106.10", "SPAC1A6.04c", "SPBC365.20c", "SPAC890.03", "SPBC8E4.02c",
        "SPAC16.01", "SPBC646.13", "SPCC4B3.12", "SPBC21C3.18", "SPAC20H4.03c",
        "SPBC17G9.09", "SPACUNK4.16c", "SPCC1020.08", "SPCC188.08c", "SPAC1783.02c",
        "SPAC13F5.04c", "SPAC23D3.03c", "SPAC8E11.05c", "SPBC1198.07c", "SPBC14C8.15",
        "SPBC18H10.18c", "SPBC1921.04c", "SPBC30D10.09c", "SPBC4B4.12c", "SPCC306.11",
        "SPCC320.03", "SPCC594.01", "SPCC594.02c", "SPCC794.03"],
    :FYPO0001309 => [
        "SPAC21E11.04", "SPCC16A11.08", "SPBC21C3.08c", "SPBC1D7.03", "SPCC1753.02c",
        "SPBC16D10.08c", "SPBC3H7.03c", "SPBC16E9.13", "SPBP4H10.11c", "SPAC17C9.02c",
        "SPAC821.07c", "SPBP35G2.11c", "SPAC806.07", "SPCC188.02", "SPCC16C4.11",
        "SPBC725.11c", "SPAC23C11.08", "SPBC3B8.02", "SPBC11B10.10c", "SPBC106.10",
        "SPAC26F1.10c", "SPBC1198.11c", "SPAC1B9.02c", "SPAC22E12.14c", "SPAC16E8.01",
        "SPCC162.12", "SPBP23A10.16", "SPBC30D10.10c", "SPAC1399.04c", "SPAC3A12.09c",
        "SPBC16E9.14c", "SPAC323.03c", "SPAC3H1.08c", "SPBP4H10.16c", "SPRRNA.47"])

end # module
