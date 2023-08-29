using DrWatson
@quickactivate 
using Revise
using Random
using SPM
using StatsBase


include(scriptsdir("cfg.jl"))

index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
h_up = 100.0
ncond = 10_000
f_tol = 1.0
x_tol = 1e-04
NOM = ARL(200)
maxrl = 10.0 * get_value(NOM)

Random.seed!(seed)
m = 500
p = 3
x = randn(m, p)

STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)


# Precompilation
bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
approximateBisectionCL(CH, nsims = 1, maxrl = 1)
saCL(CH, maxiter = 1)
time()
# end precompilation


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    fname = datadir("sims", string(typeof(STAT).name.wrapper), "bisection", "arl-$(seed)-$(nsims).jld2")
    if !isfile(fname)
        t_bisec = time()
        bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl)
        dt_bisec = time() - t_bisec
        RLs_bisec = zeros(ncond);
        for i in 1:ncond
            RLs_bisec[i] = run_sim(CH)
        end
        safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))
    end

    RLs_approx = zeros(ncond);
    for B in B_vec[i]
        fname = datadir("sims", string(typeof(STAT).name.wrapper), "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl)
            dt_approx = time() - t_approx
            for i in 1:ncond
                RLs_approx[i] = run_sim(CH)
            end
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end



#----- LLCUSUM median
NOM = QRL(200, 0.5)
STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)


# Precompilation
bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
approximateBisectionCL(CH, nsims = 1, maxrl = 1)
saCL(CH, maxiter = 1)
time()
# end precompilation


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    title = string(typeof(STAT).name.wrapper) * "-MRL"
    fname = datadir("sims", title, "bisection", "mrl-$(seed)-$(nsims).jld2")
    if !isfile(fname)
        t_bisec = time()
        bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl)
        dt_bisec = time() - t_bisec
        RLs_bisec = zeros(ncond);
        for i in 1:ncond
            RLs_bisec[i] = run_sim(CH)
        end
        safesave(fname, Dict("h" => get_h(get_limit(CH)), "mrl" => median(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))
    end

    RLs_approx = zeros(ncond);
    for B in B_vec[i]
        fname = datadir("sims", title, "approximate", "mrl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl)
            dt_approx = time() - t_approx
            for i in 1:ncond
                RLs_approx[i] = run_sim(CH)
            end
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "mrl" => median(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end