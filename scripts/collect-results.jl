using DrWatson
@quickactivate "MixedTSL"
using DataFrames
using StatsBase

results = DataFrame()
for method_ in ["bisection", "approximate"]
    METHOD = method_
    println("-------- method: $(METHOD) --------")
    title = "sims"
    for statistic in readdir(datadir(title))
        DIR = datadir(title, statistic, METHOD)
        if isdir(DIR)
            output = collect_results(DIR);
            if size(output, 1) == 0
                continue
            else
            println(describe(output, :mean, :std, :min, :median, :max))
            end
        else
            continue
        end
    end
end
# rename!(results, :ARL_mean => :AARL);
# sort!(results, [:chart, :case, :m]);
# println(results)
# println("")

# safesave(datadir("output", "ic", "results.csv"), results)
