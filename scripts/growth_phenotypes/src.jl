#=
Code for processing growth phenotype data in `process.jl`.
=#

cd(@__DIR__)

using PombeAgeingGenes, DataFrames, DataFramesMeta, JSON, CSVFiles, GZip, Statistics

const fname = ENV["POMBEAGEINGGENES"] * "/data/Jan2019_BBSRC_results"

function write_ncolonies(key::String, value::Int)
    f = ENV["POMBEAGEINGGENES"] * "/data/ncolonies.json"
    d = isfile(f) ? JSON.parsefile(f) : d = Dict{String,Int}()
    d[key] = value
    write(f, JSON.json(d))
    d
end

write_ncolonies(key::String, df::DataFrame) = write_ncolonies(key, size(df,1))

function extractcolumns()
    GZip.open(fname * ".csv.gz", "r") do input
        open(fname * "_selected_columns.csv", "w") do output
            for line = eachline(input)
                l = split(line, ',')

                replace!(l, "inf"=>"")
                replace!(l, "emty"=>"empty")

                for i = [2,3,4,5,8,9,10,13,14,15,16,22]
                    write(output, l[i], ",")
                end

                write(output, l[24], "\n")
            end
        end
    end
end

function _correct_assay_plate(c, lv, rn, old_ap, new_ap)
    df[(df[:condition] .== c) .& (df[:library_version] .== lv) .&
       (df[:repeat_number] .== rn) .& (df[:assay_plate] .== old_ap),
       :assay_plate] = new_ap
end

function correct_assay_plate(c, lv, rn, ap_A, ap_B)
    _correct_assay_plate(c, lv, rn, ap_A, "~")
    _correct_assay_plate(c, lv, rn, ap_B, ap_A)
    _correct_assay_plate(c, lv, rn, "~", ap_B)
end
