#=
Process growth phenotype data for `Oct2019_BBSRC_results.csv`
=#

using PombeAgeingGenes, DataFrames, DataFramesMeta, JSON, CSVFiles, GZip, Statistics, Plots

plotlyjs()

include("src.jl")

fname = joinpath(ENV["POMBEAGEINGGENES"], "data", "Oct2019_BBSRC_results")

#=
Extract relevent columns from the pyphe-slim output file

Columns:
    2 Column
    3 Row
    4 colony_size
    5 Colony_circularity
    8 reference_surface
    9 Colony_size_corr
    11 Scan_date
    14 condition
    15 library_version
    16 repeat_number
    17 comments
    23 ID
    25 Assay_plate
=#

GZip.open(fname * ".csv.gz", "r") do input
    open(fname * "_selected_columns.csv", "w") do output
        for line = eachline(input)
            l = split(line, ',')

            replace!(l, "inf"=>"")
            replace!(l, "emty"=>"empty")

            for i = [2,3,4,5,8,9,11,14,15,16,17,23]
                write(output, l[i], ",")
            end

            write(output, l[25], "\n")
        end
    end
end


#=
Cleaning
=#

df = DataFrame(load(fname * "_selected_columns.csv"))

rm(fname * "_selected_columns.csv")

write_ncolonies("All", df)

# Tidying up
begin
    rename!(df, [i=>Symbol(lowercase(string(i))) for i = names(df)])
    rename!(df, :colony_size_corr=>:size)

    df[df[!,:id] .== "bioneer_wt", :id] .= "wt"

    df[df[!,:comments] .== "OK", :comments] = "ok"
    df[df[!,:comments] .== "some_missing_bits", :comments] = "missing_some_bits"

    df[!,:phlox] = occursin.("phlox", df[!,:condition])
    df[!,:condition] = map(x->replace(x, "_phlox"=>""), df[!,:condition])

    for (old, new) = [
        ("YES_0.5mM_H2O2", "YES_H2O2_0.5mM"),
        ("YES_0.5mM_NaOrthovanadate", "YES_NaOrthovanadate_0.5mM"),
        ("YES_100Gr", "YES_100g"),
        ("YES_1mM_H2O2", "YES_H2O2_1mM"),
        ("YES_2.5_percent_formamide", "YES_formamide_2.5percent"),
        ("YES_Bleomycin_600ug_ml", "YES_Bleomycin_600ugml"),
        ("YES_KCl_0.5M_MMS_0.0075", "YES_KCl_0.5M_MMS_0.0075percent"),
        ("YES_MMS_0.005_percent", "YES_MMS_0.005percent"),
        ("YES_MMS_0.0075", "YES_MMS_0.0075percent"),
        ("YES_NaCl_100mM_MMS_0.075", "YES_NaCl_100mM_MMS_0.075percent"),
        ("YES_SDS_0.01_percent", "YES_SDS_0.01percent"),
        ("YES_SDS_0.02_percent", "YES_SDS_0.02percent"),
        ("YES_SDS_0.04_percent", "YES_SDS_0.04percent"),
        ("YES_Xilose_2_percent_0.1_glucose", "YES_Xilose_2percent_glucose_0.1percent"),
        ("YES_Xilose_2_percent_0.1_glucose_day2", "YES_Xilose_2percent_glucose_0.1percent_day2"),
        ("YES_benzamidine_10", "YES_benzamidine_10mM"),
        ("YES_calcofluor_2ug_ml", "YES_calcofluor_2ugml"),
        ("YES_calcofluor_2ug_ml_SDS_0.04percent", "YES_calcofluor_2ugml_SDS_0.04percent"),
        ("YES_ethanol_1percent_noglucose", "YES_ethanol_1percent_no_glucose"),
        ("YES_ethanol_5_percent", "YES_ethanol_5percent"),
        ("YES_fructose_2_percent", "YES_fructose_2percent"),
        ("YES_galactose_2_percent_0.1_glucose", "YES_galactose_2percent_glucose_0.1percent"),
        ("YES_galactose_2_percent_0.1_glucose_day2", "YES_galactose_2percent_glucose_0.1percent_day2"),
        ("YES_glycerol_MMS_0.0075", "YES_glycerol_MMS_0.0075percent"),
        ("YES_glycerol_galactose_2percent_0.1_glucose", "YES_glycerol_galactose_2percent_glucose_0.1percent"),
        ("YES_maltose_2_percent", "YES_maltose_2percent"),
        ("YES_mannitol_2_percent_01_glucose", "YES_mannitol_2percent_glucose_0.1percent"),
        ("YES_pre-rapa_rapa_100Gr", "YES_pre-rapa_rapa_100g"),
        ("YES_rapa_100Gr", "YES_rapa_100g"),
        ("YES_sucrose_2_percent", "YES_sucrose_2percent"),
        ("YES_tea_tree_0.25", "YES_tea_tree_0.25ulml"),
        ("YES_tea_tree_0.5", "YES_tea_tree_0.5ulml"),
        ]
        df[df[!,:condition] .== old, :condition] = new
    end

    df[!,:condition] = map(x->replace(x, "_percent"=>"percent"), df[!,:condition])

    dropmissing!(df, :size)

    df[!,:size] = round.(df[!,:size], digits=3)
end

# Remove rows
begin
    df = @where(df, :id .!= "grid", :id .!= "empty")

    write_ncolonies("Remove 'grid' and 'empty'", df)

    dropmissing!(df, [:size, :colony_size, :reference_surface, :colony_circularity],
                 disallowmissing=true)

    df = @where(df, :colony_size .≥ 0,
                    :size .≥ 0,
                    :reference_surface .≥ 0,
                    #.9 .< :colony_circularity .< 1.5
                    )

    write_ncolonies("Remove poor-quality colonies", df)
end

dropmissing!(df)

df = @select(df, :condition, :id, :size, :colony_size, :column, :row, :assay_plate,
                 :scan_date, :repeat_number, :library_version, :reference_surface,
                 :colony_circularity, :comments, :phlox)

save(fname * "_clean.csv", df)


#=
Normalise colony sizes
=#

df = DataFrame(load(fname * "_clean.csv"))

# For each strain, normalise the treatments to the control media (YES_32 or EMM_32)

allowmissing!(df, :size)

function controlmeans(df, media)
    x = by(@where(df, :condition .== media), :id, size = :size => mean)
    dropmissing!(x, :size)
    Dict(zip(x[:,:id], x[:,:size]))
end

const CONTROL_YES = controlmeans(df, "YES_32")
const CONTROL_EMM = controlmeans(df, "EMM_32")

for r = eachrow(df)
    d = startswith(r.condition, "YES") ? CONTROL_YES : CONTROL_EMM
    r.size = haskey(d, r.id) ? r.size / d[r.id] : missing
end

dropmissing!(df, :size)

df.size = log2.(df.size)

save(fname * "_normalised.csv", df)

#=
Process colony sizes
=#

df = DataFrame(load(fname * "_normalised.csv"))

# Select the colony size that is most different to the wild-type colony size (either
# larger or smaller)

df = by(df, [:id, :condition]) do g
    xs = g.size
    x = xs[argmax(abs.(xs))]
    (size = round(x, digits=3),)
end

#=
Wideform format for ML
=#

df = unstack(df, :id, :condition, :size)

# Impute `missing` values with mean size per condition

for col in names(df[:,Not(:id)])
    df[!,col] = coalesce.(df[!,col], mean(skipmissing(df[!,col])))
end

save(fname * "_wideform.csv", df)
