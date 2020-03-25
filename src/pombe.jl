function protein_coding_genes()
    ids = Set{String}()
    open(GzipDecompressorStream, joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
                                          "v4_2_0", "pombase_pombe_genome", "peptide.fa.gz")
    ) do io
        for line in eachline(io)
            if startswith(line, ">")
                push!(ids, split(line, ":pep")[1][2:end])
            end
        end
    end
    return ids
end
