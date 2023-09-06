using DrWatson
@quickactivate 
using Revise
using Distributions
using LinearAlgebra
using Random
using SPM
using StatsBase



include(scriptsdir("cfg.jl"))

function simulate_sensitivity(CH; target=mean, title="sims-sensitivity", statname="", h_up, ncond, maxrl, seed)
    # Precompilation
    approximateBisectionCL(CH, nsims = 1, maxrl = 1)
    time()
    # end precompilation

    if length(statname) == 0
        statname = string(typeof(STAT).name.wrapper)
    end

    Random.seed!(seed)
    f_tol_vec = [2.0, 1.0, 0.5, 0.1]
    x_tol_vec = [1e-03, 1e-05, 1e-07, 1e-09]
    B = 10000
    nsims = 10000
    for f_tol in f_tol_vec, x_tol in x_tol_vec
        RLs_approx = zeros(ncond);
        fname = datadir(title, statname, "approximate", "sim-$(seed)-$(f_tol)-$(x_tol).jld2")
        if !isfile(fname)
            t_approx = time()
            approximateBisectionCL!(CH, maxiter = 1000, nsims=nsims, B = B, maxrl = maxrl, x_tol=x_tol, f_tol=f_tol)
            dt_approx = time() - t_approx
            for i in 1:ncond
                RLs_approx[i] = run_sim(CH)
            end
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => target(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx, "x_tol" => x_tol, "f_tol" => f_tol))
        end
    end
end


#############################################################
#                          LLCUSUM                          #
#############################################################

index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
h_up = 100.0
ncond = 100_000
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

simulate_sensitivity(CH, target=mean, h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)



#----- LLCUSUM median
NOM = QRL(200, 0.5)
STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_sensitivity(CH, target=median, statname = "LLCUSUM-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)


#############################################################
#                       MEWMA 3 dim                         #
#############################################################
p = 3
seed = seeds[2] + index_sim
STAT = DiagMEWMA(Λ = fill(0.2, p))
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2Distribution(MvNormal(zeros(p), ones(p)))
NOM = ARL(200)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_sensitivity(CH, target=mean, h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)



#----- MEWMA median
NOM = QRL(200, 0.5)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_sensitivity(CH, target=median, statname="DiagMEWMA-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)



#############################################################
#                       MCUSUM 5 dim                        #
#############################################################
p = 5
seed = seeds[3] + index_sim
STAT = MCUSUM(0.25, p, 0.0, zeros(p))
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2Distribution(MvNormal(zeros(p), ones(p)))
NOM = ARL(200)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_sensitivity(CH, target=mean, statname="MCUSUM", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)



#----- AMCUSUM median
NOM = QRL(200, 0.5)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_sensitivity(CH, target=median, statname="MCUSUM-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed)
