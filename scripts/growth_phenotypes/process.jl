#=
Process growth phenotype data.
=#

include("src.jl")

#######################################
# Extract relevent columns from the pyphe-slim output file
#######################################

# Columns:
# 2 Column
# 3 Row
# 4 colony_size
# 5 Colony_circularity
# 8 reference_surface
# 9 Colony_size_corr
# 10 Scan_date
# 13 condition
# 14 library_version
# 15 repeat_number
# 16 comments
# 22 ID
# 24 Assay_plate

extractcolumns()

#######################################
# Cleaning
#######################################

df = DataFrame(load(fname * "_selected_columns.csv"))

rm(fname * "_selected_columns.csv")

write_ncolonies("All", df)

# Tidying up
begin
    rename!(df, [i=>Symbol(lowercase(string(i))) for i = names(df)])
    rename!(df, :colony_size_corr=>:size)

    df[df[:id] .== "bioneer_wt", :id] = "wt"
    df[df[:id] .== "wildtype", :id] = "wt"

    df[df[:comments] .== "OK", :comments] = "ok"
    df[df[:comments] .== "some_missing_bits", :comments] = "missing_some_bits"

    df[:phlox] = occursin.("phlox", df[:condition])

    df[:condition] = map(x->replace(x, "_phlox"=>""), df[:condition])


    for (old, new) = [
        ("YES_0.005_MMS", "YES_MMS_0.005percent"),
        ("YES_HU_10", "YES_HU_10mM"),
        ("YES_TSA_500", "YES_TSA_500nM"),
        ("EMM", "EMM_32"),
        ("YES_tunicamycin_2", "YES_tunicamycin_2ul_in_100ml"),
        ("YES_tunicamycin_1", "YES_tunicamycin_1ul_in_100ml"),
        ("YES_tunicamycin_1ulin100ml", "YES_tunicamycin_1ul_in_100ml"),
        ("YES_SDS_0.02percent", "YES_SDS_0.02_percent"),
        ("YES_SDS_0.04percent", "YES_SDS_0.04_percent"),
        ("YES_calcofluor_2ug", "YES_calcofluor_2ug_ml"),
        ("YES_calcofluor_2ug_SDS_0.04percent", "YES_calcofluor_2ug_ml_SDS_0.04percent"),
        ("YES_MMS_0.0075", "YES_MMS_0.0075percent"),
        ("YES_NaCl_100mM_MMS_0.075", "YES_NaCl_100mM_MMS_0.075percent"),
        ("YES_NaCl_100mM_SDS0.04percent", "YES_NaCl_100mM_SDS_0.04percent"),
        ("YES_LiCl_4mM_MMS_0.0075", "YES_LiCl_4mM_MMS_0.0075percent"),
        ("YES_0.1_percent_glucose", "YES_Glucose_0.1percent"),
        ("YES_0.1_percent_glucose_10days", "YES_Glucose_0.1percent_10days"),
        ("YES_0.5_percent_glucose", "YES_Glucose_0.5percent"),
        ("YES_0.5_percent_glucose_10days", "YES_Glucose_0.5percent_10days"),
        ("YES_1_percent_glucose", "YES_Glucose_1percent"),
        ("YES_1_percent_glucose_10days", "YES_Glucose_1percent_10days"),
        ("YES_3_percent_glucose_10days", "YES_Glucose_3percent_10days"),
        ("YES_0.5mM_H2O2", "YES_H2O2_0.5mM"),
        ("YES_1mM_H2O2", "YES_H2O2_1mM"),
        ("YES_0.5mM_NaOrthovanadate", "YES_NaOrthovanadate_0.5mM"),
        ("YES_2.5_percent_formamide", "YES_formamide_2.5percent"),
        ("YES_KCl_0.5M_MMS_0.0075", "YES_KCl_0.5M_MMS_0.0075percent"),
        ("YES_Galactose_0.1_glucose", "YES_Galactose_0.1percent_glucose"),
        ("YES_Glycerol_Galactose_0.1_glucose", "YES_Glycerol_Galactose_0.1percent_glucose"),
        ("YES_Glycerol_MMS_0.0075", "YES_Glycerol_MMS_0.0075percent"),
        ("YES_Glycerol_galactose_2percent_0.1_glucose", "YES_Glycerol_galactose_2percent_0.1percent_glucose"),
        ("YES_benzamidine_10", "YES_benzamidine_10mM"),
        ("YES_glucose_0.1percent", "YES_Glucose_0.1percent"),
        ("YES_tea_tree_0.25", "YES_tea_tree_0.25ul_per_ml"),
        ("YES_tea_tree_0.5", "YES_tea_tree_0.5ul_per_ml"),
        ]
        df[df[:condition] .== old, :condition] = new
    end

    df[:condition] = map(x->replace(x, "_percent"=>"percent"), df[:condition])

    df[(df[:condition] .== "YES_Xilose_2_percent_0.1_glucose") .&
       (ismissing.(df[:repeat_number])), :repeat_number] = 1.

    dropmissing!(df, :size)
    df[:size] = round.(df[:size], digits=3)

    # Rename incorrectly named plates
    correct_assay_plate("YES_MgCl2_200mM_SDS_0.4percent", "V.5.2.1", 2.2, "A", "B")
    correct_assay_plate("YES_Glycerol",	"V.5.2.2", 2.2, "A", "B")
    correct_assay_plate("YES_NaCl_100mM", "V.5.2.2", 2.4, "B", "C")
    correct_assay_plate("YES_Glucose_3percent_25C_1week", "V.5.2.2", 2.1, "A", "C")
end

# Remove rows
begin
    df = @where(df, :id .!= "grid", :id .!= "empty")

    write_ncolonies("Remove grid and empty", df)

    dropmissing!(df, [:size, :colony_size, :reference_surface, :colony_circularity],
                 disallowmissing=true)

    df = @where(df, :colony_size .> 0,
                    :size .> 0,
                    :reference_surface .> 10,
                    .9 .< :colony_circularity .< 1.5)

    write_ncolonies("Remove bad colonies", df)
end

dropmissing!(df, disallowmissing=true)

df = @select(df, :condition, :id, :size, :colony_size, :column, :row, :assay_plate,
                 :scan_date, :repeat_number, :library_version, :reference_surface,
                 :colony_circularity, :comments, :phlox)

save(fname * "_clean.csv", df)

#######################################
# Remove outliers
#######################################

df = load(GrowthPhenotypes)

removeoutliers!(df)

# These plates have more colonies than they should have
df = df[.!((df[:condition] .== "YES_formamide_2.5percent") .&
           (df[:scan_date] .== 070618)),:]

write_ncolonies("Remove outliers", df)

# Save
save(fname * "_no_outliers.csv", df)

#######################################
# Wideform format for ML
#######################################

df = load(GrowthPhenotypesNoOutliers)

# Mean sizes. NB sizes for strain-condition pairs with ≤ `nrepeats` repeats are set to `missing`
df = meansizes(df; nrepeats=1)

# For each strain, normalise the treatments to the control media (YES_32 or EMM_32)
function controlmeans(df, media)
    x = by(@where(df, :condition .== media), :id, size = :size => mean)
    dropmissing!(x, :size)
    Dict(zip(x[:id], x[:size]))
end

const CONTROL_YES = controlmeans(df, "YES_32")
const CONTROL_EMM = controlmeans(df, "EMM_32")

for r = eachrow(df)
    d = startswith(r.condition, "YES") ? CONTROL_YES : CONTROL_EMM
    r.size = haskey(d, r.id) ? r.size / d[r.id] : missing
end

# Impute `missing` values with mean size per condition
impute!(df)

df[:size] = log2.(df[:size])
df[:size] = round.(df[:size], digits=3)

# Wideform
df = unstack(df, :id, :condition, :size)
deletecols!(df, :YES_32)
deletecols!(df, :EMM_32)

# Remove low variance conditions
vars = map(i->var(skipmissing(df[i])), names(df)[2:end])
histogram(vars)
threshold = quantile(vars, 0.2)
keep_conditions = names(df)[2:end][vars .> threshold]
df = df[[:id; keep_conditions]]

# Remove low variance strains
vars = map(i->var(skipmissing(vec(Matrix(df[df[:id] .== i, 2:end])))), df[:id])
histogram(vars)
threshold = quantile(vars, 0.2)
keep_ids = df[:id][vars .> threshold]
df = @in(df, :id, keep_ids)

# Coalesce missings
impute_value = floor(minimum(skipmissing(vec(Matrix(df[2:end])))))

for col = names(df)[2:end]
    df[col] = coalesce.(df[col], impute_value)
end

save(fname * "_no_outliers_wideform.csv", df)
