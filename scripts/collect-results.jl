using DrWatson
@quickactivate 
using DataFrames
using StatsBase

results = DataFrame()
for method_ in ["bisection", "approximate"]
# for method_ in ["bisection", "approximate", "saCL"]
    METHOD = method_
    println("-------- method: $(METHOD) --------")
    title = "sims"
    for statistic in readdir(datadir(title))
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            if size(output, 1) == 0
                continue
            else
            grouped = groupby(output, [:nsims, :B, :method, :statistic])
            selector = Not(:nsims, :B, :method, :time, :statistic)
            summ = combine(grouped, :time => mean, selector .=> mean, selector .=> std, selector .=> median)
            results = vcat(results, summ)
            end
        else
            continue
        end
    end
end
sort!(results, [:statistic, :nsims, :method, :B]);
println(results)

# safesave(datadir("output", "ic", "results.csv"), results)
