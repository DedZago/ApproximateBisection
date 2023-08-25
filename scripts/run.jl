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

# Random.seed!(seed)
# m = 500
# p = 3
# x = randn(m, p)

# STAT = LLCUSUM(1.0, x)
# LIM = OneSidedFixedLimit(5.0, true)
# BOOT = Phase2(MultinomialBootstrap(STAT), x)
# CH = ControlChart(STAT, LIM, NOM, BOOT)


# # Precompilation
# bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
# approximateBisectionCL(CH, nsims = 1, maxrl = 1)
# saCL(CH, maxiter = 1)
# time()
# # end precompilation


# # Random.seed!(seed)
# # t_sacl = time()
# # saCL!(CH, gamma = 0.01, maxiter = 1000000, verbose=false)
# # dt_sacl = time() - t_sacl
# # RLs_sacl = zeros(ncond);
# # for i in 1:ncond
# #     RLs_sacl[i] = run_sim(CH)
# # end
# # safesave(datadir("sims", string(typeof(STAT).name.wrapper), "saCL", "arl-$(seed).jld2"), Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_sacl), "nsims" => 0, "B" => 0, "time" => dt_sacl))



# Random.seed!(seed)
# nsims_vec = [1000, 10000]
# B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
# for i in 1:length(nsims_vec)
#     nsims = nsims_vec[i] 
#     fname = datadir("sims", string(typeof(STAT).name.wrapper), "bisection", "arl-$(seed)-$(nsims).jld2")
#     if !isfile(fname)
#         t_bisec = time()
#         bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl)
#         dt_bisec = time() - t_bisec
#         RLs_bisec = zeros(ncond);
#         for i in 1:ncond
#             RLs_bisec[i] = run_sim(CH)
#         end
#         safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))
#     end

#     RLs_approx = zeros(ncond);
#     for B in B_vec[i]
#         fname = datadir("sims", string(typeof(STAT).name.wrapper), "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
#         if !isfile(fname)
#             t_approx = time()
#             approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl)
#             dt_approx = time() - t_approx
#             for i in 1:ncond
#                 RLs_approx[i] = run_sim(CH)
#             end
#             safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => mean(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
#         end
#     end
# end

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


# #--- MULTIPLE CHART WITH LLCUSUM
# using Distributions
# seed = seeds[3] + index_sim
# Random.seed!(seed)
# statistic_name = "MultipleLLCUSUM"
# m = 500
# p = 2
# x = randn(m, p)
# NM = ARL(200)
# STAT1 = LLCUSUM(0.1, x)
# STAT2 = LLCUSUM(0.5, x)
# LIM1 = OneSidedFixedLimit(1.0, true)
# LIM2 = OneSidedFixedLimit(1.0, true)
# PH2 = Phase2(MultinomialBootstrap(STAT1), x)
# CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NM, PH2)

# # # Precompilation
# approximateBisectionCL(CH, nsims = 5, maxrl = 5, f_tol = f_tol, x_tol = x_tol)
# saCL(CH, maxiter=1)
# time()
# # # end precompilation


# Random.seed!(seed)
# fname = datadir("sims", statistic_name, "saCL", "arl-$(seed).jld2")
# if !isfile(fname)
#     t_sacl = time()
#     saCL!(CH, Amax=2.5, gamma = 0.01, maxiter = 10000000, verbose=false)
#     dt_sacl = time() - t_sacl
#     # RLs_sacl = zeros(ncond);
#     tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
#     ARL0 = mean(minimum(x[:rl]) for x in tmp)
#     ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
#     safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
# end


# Random.seed!(seed)
# nsims_vec = [1000, 10000]
# B_vec = [[nsims_vec[1]], [nsims_vec[2]]]
# for i in 1:length(nsims_vec)
#     nsims = nsims_vec[i] 
#     for B in B_vec[i]
#         fname = datadir("sims", statistic_name, "approximate", "arl-$(seed)-$(nsims)-$(B).jld2")
#         if !isfile(fname)
#             t_approx = time()
#             approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
#             dt_approx = time() - t_approx
#             tmp = [run_sim_sa(CH, maxiter=10^4, delta=0.0) for _ in 1:ncond]
#             ARL0 = mean(minimum(x[:rl]) for x in tmp)
#             ARL_js = [mean(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
#             safesave(fname, Dict("h" => get_h(get_limit(CH)), "arl" => ARL0, "arl_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
#         end
#     end
# end




#--- MULTIPLE CHART WITH LLCUSUM SMALLER k
using Distributions
seed = seeds[3] + index_sim
Random.seed!(seed)
statistic_name = "MultipleLLCUSUMsmall"
m = 500
p = 2
x = randn(m, p)
NM = ARL(200)
STAT1 = LLCUSUM(0.1, x)
STAT2 = LLCUSUM(0.01, x)
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
