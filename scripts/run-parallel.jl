using DrWatson
@quickactivate 
using Revise
using Distributions
using LinearAlgebra
using Random
using SPM
using StatsBase



include(scriptsdir("cfg.jl"))
index_sim = parsed_args["index"]

function simulate_parallel(CH; target=mean, title="sims-parallel", statname="", h_up, ncond, maxrl, seed, f_tol, x_tol, parallel=true)
    # Precompilation
    approximateBisectionCL(CH, nsims = 2, maxrl = 3, parallel=true)
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
        approximateBisectionCL!(CH, maxiter = 1000, nsims=nsims, B = B, maxrl = maxrl, x_tol=x_tol, f_tol=f_tol, parallel=parallel)
        dt_approx = time() - t_approx
        tmp = Vector{NamedTuple{(:rl, :rlPlus, :rlMinus), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}(undef, ncond)
        Threads.@threads for i in 1:ncond
            tmp[i] = run_sim_sa(CH, maxiter=10^5, delta=0.0)
        end
        ARL0 = target(minimum(x[:rl]) for x in tmp)
        ARL_js = [target(x[:rl][l] for x in tmp) for l in 1:length(get_statistic(CH))]
        safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => ARL0, "target_j" => ARL_js, "nsims" => nsims, "B" => B, "time" => dt_approx))
    end
end


seed = seeds[1] + index_sim
h_up = 100.0
ncond = 10_000
f_tol = 1.0
x_tol = 1e-06
NOM = ARL(200)
maxrl = 10.0 * get_value(NOM)


############################################################
#                      MULTIPLE RSADA                      #
############################################################
println("Multiple RSADA")
using Distributions
seed = seeds[1] + index_sim
statname = "RSADA"
NOM = ARL(200)
p = 200
q = 20
x = randn(1,p)
STAT1 = RSADA(0.5, 1.5, q, x)
STAT2 = RSADA(0.75, 2.0, q, x)
STAT3 = RSADA(1.0, 2.5, q, x)
STAT4 = RSADA(1.25, 3.0, q, x)
LIM1 = OneSidedFixedLimit(1.0, true)
LIM2 = OneSidedFixedLimit(1.0, true)
LIM3 = OneSidedFixedLimit(1.0, true)
LIM4 = OneSidedFixedLimit(1.0, true)
PH2 = Phase2Distribution(MvNormal(zeros(p), ones(p)))
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)

simulate_parallel(CH, target=mean, h_up=h_up, statname = statname, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol, title="sims-parallel-8cores")
GC.gc(true)


println("Multiple RSADA median optimization")
NOM = QRL(200, 0.5)
# maxrl = 20.0 * get_value(NOM)
CH = ControlChart([STAT1, STAT2, STAT3, STAT4], [LIM1, LIM2, LIM3, LIM4], NOM, PH2)

simulate_parallel(CH, target=median, statname = statname * "-median", h_up=h_up, ncond=ncond, maxrl=maxrl, seed=seed, f_tol=f_tol, x_tol=x_tol, title="sims-parallel-8cores")
