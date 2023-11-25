using DrWatson
@quickactivate 
using DataFrames
using StatsBase
using CSVFiles
using Latexify

title = "sims"
for statistic in readdir(datadir(title))
    results = DataFrame()
    for method_ in ["bisection", "approximate"]
        METHOD = method_
        println("-------- method: $(METHOD) --------")
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            # output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            if :target_j in propertynames(output)
                expanded_targetj = Matrix(hcat(output.target_j...)')
                df_targetj = DataFrame(expanded_targetj, Symbol.("target_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                expanded_h = Matrix(hcat(output.h...)')
                df_h = DataFrame(expanded_h, Symbol.("h_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                output = output[:, Not(:target_j, :h)]
                output = hcat(output, df_targetj, df_h)
            end
            results = vcat(results, output)
        else
            continue
        end
    end
    grouped = groupby(results, [:nsims, :method, :B])
    selector = Not(:nsims, :method, :B)
    summ = combine(grouped, selector .=> median, selector .=> iqr)
    num_columns = [i for i in names(summ) if Base.nonmissingtype(eltype(summ[!,i])) <: Number]
    summ_numerical = summ[!, num_columns]
    summ_non_numerical = summ[!, Not(num_columns)]
    summ_numerical = round.(summ_numerical, digits=3)
    summ_final = hcat(summ_non_numerical, summ_numerical)
    sort!(summ_final, [:nsims,:method,:B])
    summ_final_transpose = DataFrame([[names(summ_final)]; collect.(eachrow(summ_final))], [:column; Symbol.(axes(summ_final, 1))])
    println(summ_final)
    # safesave(datadir("output", title, statistic, "results.csv"), summ)
    # @show latexify(summ_final)
    safesave(datadir("output", title, statistic, "results_rounded.csv"), summ_final)
    write(datadir("output", title, statistic, "results.tex"), latexify(summ_final, env=:table))
end


title = "sims-sacl"
for statistic in readdir(datadir(title))
    results = DataFrame()
    for method_ in ["approximate", "saCL"]
        METHOD = method_
        println("-------- method: $(METHOD) --------")
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            # output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            if :target_j in propertynames(output)
                expanded_targetj = Matrix(hcat(output.target_j...)')
                df_targetj = DataFrame(expanded_targetj, Symbol.("target_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                expanded_h = Matrix(hcat(output.h...)')
                df_h = DataFrame(expanded_h, Symbol.("h_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                output = output[:, Not(:target_j, :h)]
                output = hcat(output, df_targetj, df_h)
            end
            results = vcat(results, output)
        else
            continue
        end
    end
    grouped = groupby(results, [:nsims, :method, :B])
    selector = Not(:nsims, :method, :B)
    summ = combine(grouped, selector .=> median, selector .=> iqr)
    num_columns = [i for i in names(summ) if Base.nonmissingtype(eltype(summ[!,i])) <: Number]
    summ_numerical = summ[!, num_columns]
    summ_non_numerical = summ[!, Not(num_columns)]
    summ_numerical = round.(summ_numerical, digits=4)
    summ_final = hcat(summ_non_numerical, summ_numerical)
    sort!(summ_final, [:nsims,:method,:B])
    summ_final = summ_final[!, Not(:nsims, :B)]
    summ_final_transpose = DataFrame([[names(summ_final)]; collect.(eachrow(summ_final))], [:column; Symbol.(axes(summ_final, 1))])
    summ_final_transpose = summ_final_transpose[2:nrow(summ_final_transpose), 1:ncol(summ_final_transpose)]
    # col_names = names(summ_final)
    # col_names = [replace(x, "median" => "") for x in col_names]
    # col_names = [replace(x, "iqr" => "") for x in col_names]
    rename!(summ_final_transpose, Symbol.([" ", "SA", "BA-Bisection"]))
    println(summ_final_transpose)
    # safesave(datadir("output", title, statistic, "results.csv"), summ_final_transpose)
    # @show latexify(summ_final_transpose)
    safesave(datadir("output", title, statistic, "results_rounded.csv"), summ_final)
    write(datadir("output", title, statistic, "results.tex"), latexify(summ_final_transpose, env=:table))
end



title = "sims-parallel"
for statistic in readdir(datadir(title))
    results = DataFrame()
    for method_ in ["approximate"]
        METHOD = method_
        println("-------- method: $(METHOD) --------")
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            # output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            if :target_j in propertynames(output)
                expanded_targetj = Matrix(hcat(output.target_j...)')
                df_targetj = DataFrame(expanded_targetj, Symbol.("target_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                expanded_h = Matrix(hcat(output.h...)')
                df_h = DataFrame(expanded_h, Symbol.("h_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                output = output[:, Not(:target_j, :h)]
                output = hcat(output, df_targetj, df_h)
            end
            results = vcat(results, output)
        else
            continue
        end
    end
    grouped = groupby(results, [:nsims, :method, :B])
    selector = Not(:nsims, :method, :B)
    summ = combine(grouped, selector .=> median, selector .=> iqr)
    num_columns = [i for i in names(summ) if Base.nonmissingtype(eltype(summ[!,i])) <: Number]
    summ_numerical = summ[!, num_columns]
    summ_non_numerical = summ[!, Not(num_columns)]
    summ_numerical = round.(summ_numerical, digits=4)
    summ_final = hcat(summ_non_numerical, summ_numerical)
    sort!(summ_final, [:nsims,:method,:B])
    summ_final = summ_final[!, Not(:nsims, :B)]
    summ_final_transpose = DataFrame([[names(summ_final)]; collect.(eachrow(summ_final))], [:column; Symbol.(axes(summ_final, 1))])
    summ_final_transpose = summ_final_transpose[2:nrow(summ_final_transpose), 1:ncol(summ_final_transpose)]
    # col_names = names(summ_final)
    # col_names = [replace(x, "median" => "") for x in col_names]
    # col_names = [replace(x, "iqr" => "") for x in col_names]
    rename!(summ_final_transpose, Symbol.([" ", "BA-Bisection"]))
    println(summ_final_transpose)
    # safesave(datadir("output", title, statistic, "results.csv"), summ_final_transpose)
    # @show latexify(summ_final_transpose)
    safesave(datadir("output", title, statistic, "results_rounded.csv"), summ_final)
    write(datadir("output", title, statistic, "results.tex"), latexify(summ_final_transpose, env=:table))
end


title = "sims-nonparallel"
for statistic in readdir(datadir(title))
    results = DataFrame()
    for method_ in ["approximate"]
        METHOD = method_
        println("-------- method: $(METHOD) --------")
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            # output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            if :target_j in propertynames(output)
                expanded_targetj = Matrix(hcat(output.target_j...)')
                df_targetj = DataFrame(expanded_targetj, Symbol.("target_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                expanded_h = Matrix(hcat(output.h...)')
                df_h = DataFrame(expanded_h, Symbol.("h_" .* [string(i) for i in 1:size(expanded_targetj, 2)]))
                output = output[:, Not(:target_j, :h)]
                output = hcat(output, df_targetj, df_h)
            end
            results = vcat(results, output)
        else
            continue
        end
    end
    grouped = groupby(results, [:nsims, :method, :B])
    selector = Not(:nsims, :method, :B)
    summ = combine(grouped, selector .=> median, selector .=> iqr)
    num_columns = [i for i in names(summ) if Base.nonmissingtype(eltype(summ[!,i])) <: Number]
    summ_numerical = summ[!, num_columns]
    summ_non_numerical = summ[!, Not(num_columns)]
    summ_numerical = round.(summ_numerical, digits=4)
    summ_final = hcat(summ_non_numerical, summ_numerical)
    sort!(summ_final, [:nsims,:method,:B])
    summ_final = summ_final[!, Not(:nsims, :B)]
    summ_final_transpose = DataFrame([[names(summ_final)]; collect.(eachrow(summ_final))], [:column; Symbol.(axes(summ_final, 1))])
    summ_final_transpose = summ_final_transpose[2:nrow(summ_final_transpose), 1:ncol(summ_final_transpose)]
    # col_names = names(summ_final)
    # col_names = [replace(x, "median" => "") for x in col_names]
    # col_names = [replace(x, "iqr" => "") for x in col_names]
    rename!(summ_final_transpose, Symbol.([" ", "BA-Bisection"]))
    println(summ_final_transpose)
    # safesave(datadir("output", title, statistic, "results.csv"), summ_final_transpose)
    # @show latexify(summ_final_transpose)
    safesave(datadir("output", title, statistic, "results_rounded.csv"), summ_final)
    write(datadir("output", title, statistic, "results.tex"), latexify(summ_final_transpose, env=:table))
end

