using DrWatson
@quickactivate 
using DataFrames
using StatsBase

title = "sims"
for statistic in readdir(datadir(title))
    results = DataFrame()
    for method_ in ["bisection", "approximate", "saCL"]
        METHOD = method_
        println("-------- method: $(METHOD) --------")
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            output.method = [METHOD for _ in eachrow(output)]
            # output.statistic = [statistic for _ in eachrow(output)]
            output = output[:, Not(:path)]
            results = vcat(results, output)
        else
            continue
        end
    end
    grouped = groupby(results, [:nsims, :method, :B])
    selector = Not(:nsims, :B, :method)
    summ = combine(grouped, selector .=> mean, selector .=> std)
    sort!(summ, [:nsims, :B, :method])
    println(summ)
end


