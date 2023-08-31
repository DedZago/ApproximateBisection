using DrWatson
@quickactivate 
using Revise
using Random
using SPM
using StatsBase


include(scriptsdir("cfg.jl"))

function simulate_control_chart_sacl(CH; target=mean, title="sims-sacl", statname="", ncond, f_tol, x_tol, maxrl, seed, gamma=0.01, Amax)

    # Precompilation
    approximateBisectionCL(CH, nsims = 1, maxrl = 1)
    saCL(CH, maxiter = 1)
    time()
    # end precompilation

    if length(statname) == 0
        statname = string(typeof(STAT).name.wrapper)
    end

    Random.seed!(seed)

    fname = datadir(title, statname, "saCL", "sim-$(seed).jld2")
    if !isfile(fname)
        t_sacl = time()
        saCL!(CH, Amax=Amax, gamma=gamma, maxiter = 1000000, verbose=false)
        dt_sacl = time() - t_sacl
        # RLs_sacl = zeros(ncond);
        tmp = [run_sim_sa(CH, maxiter=10^5, delta=0.0) for _ in 1:ncond]
        ARL0 = target(minimum(x[:rl]) for x in tmp)
        ARL_js = [target(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
        safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => ARL0, "target_j" => ARL_js, "nsims" => 0, "B" => 0, "time" => dt_sacl))
    end

    Random.seed!(seed)
    nsims_vec = [10000]
    B_vec = [10000]
    for i in 1:length(nsims_vec)
        nsims = nsims_vec[i] 
        for B in B_vec[i]
            fname = datadir(title, statname, "approximate", "sim-$(seed)-$(nsims)-$(B).jld2")
            if !isfile(fname)
                t_approx = time()
                approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, f_tol = f_tol, x_tol = x_tol)
                dt_approx = time() - t_approx
                tmp = [run_sim_sa(CH, maxiter=10^5, delta=0.0) for _ in 1:ncond]
                ARL0 = target(minimum(x[:rl]) for x in tmp)
                ARL_js = [target(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
                safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => ARL0, "target_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
            end
        end
    end

end


index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
h_up = 100.0
ncond = 10_000
f_tol = 1.0
x_tol = 1e-04
NOM = ARL(200)
maxrl = 10.0 * get_value(NOM)
gamma = 0.01

#--- MULTIPLE CHART WITH LLCUSUM
println("Multiple LLCUSUM")
using Distributions
seed = seeds[1] + index_sim
Random.seed!(seed)
statname = "MultipleLLCUSUM"
m = 500
p = 2
x = randn(m, p)
NOM = ARL(200)
STAT1 = LLCUSUM(0.1, x)
STAT2 = LLCUSUM(0.5, x)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2(MultinomialBootstrap(STAT1), x)
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NOM, PH2)

simulate_control_chart_sacl(CH, target=mean, statname=statname, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 2.5)

println("Multiple LLCUSUM median optimization")
NOM = QRL(200, 0.5)
maxrl = 20.0 * get_value(NOM)
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NOM, PH2)
simulate_control_chart_sacl(CH, target=median, statname=statname*"-median", ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 2.0)


#--- MULTIPLE CHART WITH T2-MEWMA
println("T2-MEWMA")
using Distributions
using LinearAlgebra
seed = seeds[2] + index_sim
statname = "T2-MEWMA"
NOM = ARL(200)
p = 3
maxrl = 10.0 * get_value(NOM)
STAT1 = MShewhart(μ = zeros(p), Σ_m1 = diagm(ones(p)))
STAT2 = DiagMEWMA(Λ = [0.2 for _ in 1:p])
LIM1 = OneSidedFixedLimit(2.0, true)
LIM2 = OneSidedFixedLimit(2.0, true)
PH2 = Phase2Distribution(MvNormal(zeros(p), ones(p)))
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NOM, PH2)

simulate_control_chart_sacl(CH, target=mean, statname=statname, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 5.0)

println("T2-MEWMA median optimization")
NOM = QRL(200, 0.5)
maxrl = 20.0 * get_value(NOM)
CH = ControlChart([STAT1, STAT2], [LIM1, LIM2], NOM, PH2)
simulate_control_chart_sacl(CH, target=median, statname=statname*"-median", ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 3.0)


#--- MULTIPLE CHART WITH EWMA
println("Multiple EWMA")
using Distributions
seed = seeds[3] + index_sim
statname = "MultipleEWMA"
NOM = ARL(200)
maxrl = 10.0 * get_value(NOM)
STAT1 = EWMA(λ = 0.2)
STAT2 = EWMA(λ = 0.5)
STAT3 = EWMA(λ = 0.05)
STAT4 = EWMA(λ = 0.1)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
LIM3 = OneSidedFixedLimit(1.0, true)
LIM4 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2Distribution(Normal(0,1))
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)

simulate_control_chart_sacl(CH, target=mean, statname=statname, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 1.0)


println("Multiple EWMA median optimization")
NOM = QRL(200, 0.5)
maxrl = 20.0 * get_value(NOM)
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)
simulate_control_chart_sacl(CH, target=median, statname=statname*"-median", ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed, gamma = gamma, Amax = 1.0)

