using PombeAgeingGenes, DataFrames, DataFramesMeta, JSON, CSVFiles, GZip, Statistics, Plots

function write_ncolonies(key::String, value::Int)
    f = joinpath(ENV["POMBEAGEINGGENES"], "data", "ncolonies.json")
    d = isfile(f) ? JSON.parsefile(f) : d = Dict{String,Int}()
    d[key] = value
    write(f, JSON.json(d))
    return d
end

write_ncolonies(key::String, df::DataFrame) = write_ncolonies(key, size(df,1))

function _correct_assay_plate(c, lv, rn, old_ap, new_ap)
    df[(df[!,:condition] .== c) .& (df[!,:library_version] .== lv) .&
       (df[!,:repeat_number] .== rn) .& (df[!,:assay_plate] .== old_ap),
       :assay_plate] = new_ap
end

function correct_assay_plate(c, lv, rn, ap_A, ap_B)
    _correct_assay_plate(c, lv, rn, ap_A, "~")
    _correct_assay_plate(c, lv, rn, ap_B, ap_A)
    _correct_assay_plate(c, lv, rn, "~", ap_B)
end
