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

#--- MULTIPLE CHART WITH EWMA
using Distributions
seed = seeds[2] + index_sim
statistic_name = "MultipleEWMA"
NM = ARL(200)
STAT1 = EWMA(λ = 0.2)
STAT2 = EWMA(λ = 0.5)
STAT3 = EWMA(λ = 0.05)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(2.0, true)
LIM3 = OneSidedFixedLimit(0.1, true)
PH1 = Phase2Distribution(Normal(0,1))
CH = ControlChart([STAT1, STAT2, STAT3], [LIM1, LIM2, LIM3], NM, PH1)

# # Precompilation
approximateBisectionCL(CH, nsims = 5, maxrl = 5, f_tol = f_tol, x_tol = x_tol)
saCL(CH, maxiter=1)
time()
# # end precompilation


Random.seed!(seed)
fname = datadir("sims", statistic_name, "saCL", "arl-$(seed).jld2")
if !isfile(fname)
    t_sacl = time()
    saCL!(CH, Amax=10.0, gamma = 0.01, maxiter = 10000000, verbose=false)
    dt_sacl = time() - t_sacl
    # RLs_sacl = zeros(ncond);
    tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
    ARL0 = mean(minimum(x[:rl]) for x in tmp)
    ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
    safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
end


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[nsims_vec[1]], [nsims_vec[2]]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    for B in B_vec[i]
        fname = datadir("sims", statistic_name, "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
            dt_approx = time() - t_approx
            tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
            ARL0 = mean(minimum(x[:rl]) for x in tmp)
            ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end



#--- MULTIPLE CHART WITH LLCUSUM
using Distributions
seed = seeds[3] + index_sim
Random.seed!(seed)
statistic_name = "MultipleLLCUSUM"
m = 500
p = 2
x = randn(m, p)
NM = ARL(200)
STAT1 = LLCUSUM(0.1, x)
STAT2 = LLCUSUM(0.5, x)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2(MultinomialBootstrap(STAT1), x)
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NM, PH2)

# # Precompilation
approximateBisectionCL(CH, nsims = 5, maxrl = 5, f_tol = f_tol, x_tol = x_tol)
saCL(CH, maxiter=1)
time()
# # end precompilation


Random.seed!(seed)
fname = datadir("sims", statistic_name, "saCL", "arl-$(seed).jld2")
if !isfile(fname)
    t_sacl = time()
    saCL!(CH, Amax=2.5, gamma = 0.01, maxiter = 10000000, verbose=false)
    dt_sacl = time() - t_sacl
    # RLs_sacl = zeros(ncond);
    tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
    ARL0 = mean(minimum(x[:rl]) for x in tmp)
    ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
    safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
end


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[nsims_vec[1]], [nsims_vec[2]]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    for B in B_vec[i]
        fname = datadir("sims", statistic_name, "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
            dt_approx = time() - t_approx
            tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
            ARL0 = mean(minimum(x[:rl]) for x in tmp)
            ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end



#--- MULTIPLE CHART WITH T2-MEWMA
using Distributions
using LinearAlgebra
seed = seeds[4] + index_sim
statistic_name = "T2-MEWMA"
NM = ARL(200)
p = 3
STAT1 = MShewhart(μ = zeros(p), Σ = diagm(ones(p)))
STAT2 = DiagMEWMA(λ = 0.2)
STAT2 = EWMA(λ = 0.5)
STAT3 = EWMA(λ = 0.05)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(2.0, true)
LIM3 = OneSidedFixedLimit(0.1, true)
PH1 = Phase2Distribution(Normal(0,1))
CH = ControlChart([STAT1, STAT2, STAT3], [LIM1, LIM2, LIM3], NM, PH1)

# # Precompilation
approximateBisectionCL(CH, nsims = 5, maxrl = 5, f_tol = f_tol, x_tol = x_tol)
saCL(CH, maxiter=1)
time()
# # end precompilation


Random.seed!(seed)
fname = datadir("sims", statistic_name, "saCL", "arl-$(seed).jld2")
if !isfile(fname)
    t_sacl = time()
    saCL!(CH, Amax=10.0, gamma = 0.01, maxiter = 10000000, verbose=false)
    dt_sacl = time() - t_sacl
    # RLs_sacl = zeros(ncond);
    tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
    ARL0 = mean(minimum(x[:rl]) for x in tmp)
    ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
    safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
end


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[nsims_vec[1]], [nsims_vec[2]]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    for B in B_vec[i]
        fname = datadir("sims", statistic_name, "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
            dt_approx = time() - t_approx
            tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
            ARL0 = mean(minimum(x[:rl]) for x in tmp)
            ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end



#--- MULTIPLE CHART WITH LLCUSUM
using Distributions
seed = seeds[3] + index_sim
Random.seed!(seed)
statistic_name = "MultipleLLCUSUM"
m = 500
p = 2
x = randn(m, p)
NM = ARL(200)
STAT1 = LLCUSUM(0.1, x)
STAT2 = LLCUSUM(0.5, x)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2(MultinomialBootstrap(STAT1), x)
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NM, PH2)

# # Precompilation
approximateBisectionCL(CH, nsims = 5, maxrl = 5, f_tol = f_tol, x_tol = x_tol)
saCL(CH, maxiter=1)
time()
# # end precompilation


Random.seed!(seed)
fname = datadir("sims", statistic_name, "saCL", "arl-$(seed).jld2")
if !isfile(fname)
    t_sacl = time()
    saCL!(CH, Amax=2.5, gamma = 0.01, maxiter = 10000000, verbose=false)
    dt_sacl = time() - t_sacl
    # RLs_sacl = zeros(ncond);
    tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
    ARL0 = mean(minimum(x[:rl]) for x in tmp)
    ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
    safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
end


Random.seed!(seed)
nsims_vec = [1000, 10000]
B_vec = [[nsims_vec[1]], [nsims_vec[2]]]
for i in 1:length(nsims_vec)
    nsims = nsims_vec[i] 
    for B in B_vec[i]
        fname = datadir("sims", statistic_name, "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
            dt_approx = time() - t_approx
            tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
            ARL0 = mean(minimum(x[:rl]) for x in tmp)
            ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
        end
    end
end