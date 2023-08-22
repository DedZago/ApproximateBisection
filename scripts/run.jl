using DrWatson
@quickactivate "ApproximateBisection"
using Revise
using Random
using SPM
using StatsBase


include(scriptsdir("cfg.jl"))

index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
Random.seed!(seed)
m = 500
p = 3
x = randn(m, p)

STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
NOM = ARL(200)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)
maxrl = 10.0 * get_value(NOM)
h_up = 100.0


# Precompilation
bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
approximateBisectionCL(CH, nsims = 1, maxrl = 1)
saCL(CH, maxiter = 1)
time()
# end precompilation


Random.seed!(seed)
t_sacl = time()
saCL!(CH, gamma = 0.01, maxiter = 10000000, verbose=false)
dt_sacl = time() - t_sacl
RLs_sacl = zeros(ncond);
for i in 1:ncond
    RLs_sacl[i] = run_sim(CH)
end
safesave(datadir("sims", string(typeof(STAT).name.wrapper), "saCL", "arl-$(seed)-$(nsims).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_sacl), "nsims" => nsims, "B" => nsims, "time" => dt_sacl))



Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
ncond = 100_000
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    t_bisec = time()
    bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl)
    dt_bisec = time() - t_bisec
    RLs_bisec = zeros(ncond);
    for i in 1:ncond
        RLs_bisec[i] = run_sim(CH)
    end
    safesave(datadir("sims", string(typeof(STAT).name.wrapper), "bisection", "arl-$(seed)-$(nsims).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))

    RLs_approx = zeros(ncond);
    for B in B_vec[i]
        t_approx = time()
        approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl)
        dt_approx = time() - t_approx
        for i in 1:ncond
            RLs_approx[i] = run_sim(CH)
        end
        safesave(datadir("sims", string(typeof(STAT).name.wrapper), "approximate", "arl-$(seed)-$(nsims)-$(B).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
    end
end


#--- AMCUSUM
STAT = AMCUSUM(0.2, x)
LIM = OneSidedFixedLimit(5.0, true)
NOM = ARL(200)
BOOT = Phase2(Bootstrap(), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)

# Precompilation
bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
approximateBisectionCL(CH, nsims = 1, maxrl = 1)
time()
# end precompilation

Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
ncond = 100_000
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    t_bisec = time()
    bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl)
    dt_bisec = time() - t_bisec
    RLs_bisec = zeros(ncond);
    for i in 1:ncond
        RLs_bisec[i] = run_sim(CH)
    end
    safesave(datadir("sims", string(typeof(STAT).name.wrapper), "bisection", "arl-$(seed)-$(nsims).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))

    RLs_approx = zeros(ncond);
    for B in B_vec[i]
        t_approx = time()
        approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl)
        dt_approx = time() - t_approx
        for i in 1:ncond
            RLs_approx[i] = run_sim(CH)
        end
        safesave(datadir("sims", string(typeof(STAT).name.wrapper), "approximate", "arl-$(seed)-$(nsims)-$(B).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
    end
end