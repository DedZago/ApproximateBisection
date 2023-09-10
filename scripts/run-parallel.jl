using DrWatson
@quickactivate 
using Revise
using Distributions
using LinearAlgebra
using Random
using SPM
using StatsBase



include(scriptsdir("cfg.jl"))

function simulate_parallel(CH; target=mean, title="sims-parallel", statname="", h_up, ncond, maxrl, seed, f_tol, x_tol)
    # Precompilation
    approximateBisectionCL(CH, nsims = 1, maxrl = 1, parallel=true)
    time()
    # end precompilation

    if length(statname) == 0
        statname = string(typeof(STAT).name.wrapper)
    end

    Random.seed!(seed)
    B = 10000
    nsims = 10000
    RLs_approx = zeros(ncond);
    fname = datadir(title, statname, "approximate", "sim-$(seed).jld2")
    @show fname
    if !isfile(fname)
        t_approx = time()
        approximateBisectionCL!(CH, maxiter = 1000, nsims=nsims, B = B, maxrl = maxrl, x_tol=x_tol, f_tol=f_tol, parallel=true)
        dt_approx = time() - t_approx
        Threads.@threads for i in 1:ncond
            RLs_approx[i] = run_sim(CH)
        end
        safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => target(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx, "x_tol" => x_tol, "f_tol" => f_tol))
    end
end


#############################################################
#                          LLCUSUM                          #
#############################################################

index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
h_up = 100.0
ncond = 100_000
f_tol = 1.0
x_tol = 1e-06
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

simulate_parallel(CH, target=mean, h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol)
GC.gc(true)


#----- LLCUSUM median
NOM = QRL(200, 0.5)
STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_parallel(CH, target=median, statname = "LLCUSUM-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol)


############################################################
#                      MULTIPLE RSADA                      #
############################################################
println("Multiple RSADA")
using Distributions
seed = seeds[1] + index_sim
statname = "RSADA"
NOM = ARL(200)
p = 100
q = 15
x = randn(1,p)
STAT1 = RSADA(0.15, 0.5, q, x)
STAT2 = RSADA(0.3, 1.0, q, x)
STAT3 = RSADA(0.5, 1.5, q, x)
STAT4 = RSADA(1.0, 2.0, q, x)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
LIM3 = OneSidedFixedLimit(1.0, true)
LIM4 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2Distribution(MvNormal(zeros(p), ones(p)))
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)

simulate_parallel(CH, target=mean, h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol)
GC.gc(true)


println("Multiple RSADA median optimization")
NOM = QRL(200, 0.5)
# maxrl = 20.0 * get_value(NOM)
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)

simulate_parallel(CH, target=median, statname = statname * "-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol)