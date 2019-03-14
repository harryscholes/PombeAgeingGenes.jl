#=
Process growth phenotype data.
=#

include("src.jl")

#######################################
# Extract relevent columns from the pyphe-slim output file
#######################################

# Columns:
#     2 Column
#     3 Row
#     4 colony_size
#     5 Colony_circularity
#     8 reference_surface
#     9 Colony_size_corr
#     10 Scan_date
#     13 condition
#     14 library_version
#     15 repeat_number
#     16 comments
#     22 ID
#     24 Assay_plate

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

    df[(df[:condition] .== "YES_Xilose_2_percent_0.1_glucose") .&
       (ismissing.(df[:repeat_number])), :repeat_number] = 1.

    df[:size] = round.(df[:size], digits=2)

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

    dropmissing!(df, [:size, :colony_size, :reference_surface, :colony_circularity], disallowmissing=true)

    df = @where(df, :colony_size .> 0,
                    :size .> 0,
                    :reference_surface .> 10,
                    .9 .< :colony_circularity .< 1.5)

    write_ncolonies("Remove bad colonies", df)
end

dropmissing!(df, disallowmissing=true)

df = @select(df, :condition, :id, :size, :colony_size, :column, :row, :assay_plate, :scan_date, :repeat_number, :library_version, :reference_surface, :colony_circularity, :comments, :phlox)

save(fname * "_clean.csv", df)

#######################################
# Remove outliers
#######################################

df = load(GrowthPhenotypes)

removeoutliers!(df)

# These plates have more colonies than they should have
df = df[.!((df[:condition] .== "YES_2.5_percent_formamide") .&
           (df[:scan_date] .== 070618)),:]

write_ncolonies("Remove outliers", df)

# Clamp large values
df[df[:size] .> 2, :size] = 2

# Save
save(fname * "_no_outliers.csv", df)

#######################################
# Wideform format for ML
#######################################

df = load(GrowthPhenotypesNoOutliers)

# Mean sizes. NB sizes for strain-condition pairs with ≤ `nrepeats` repeats are set to NaN
df = meansizes(df; nrepeats=2)

# Remove columns with lots of NaNs
df_nnan = by(df, :condition, nnan = :size => x->count(isnan.(x)))
conditions_to_delete = @where(df_nnan, :nnan .> 3000)[:condition]
deleterows!(df, map(x->x ∈ conditions_to_delete, df[:condition]))

# Impute NaNs with mean size per condition
impute!(df)

# Wideform
df = unstack(df, :id, :condition, :size)

# Coalesce missings
for col = names(df)[2:end]
    df[col] = coalesce.(df[col], 0.)
end

save(fname * "_no_outliers_wideform.csv", df)
