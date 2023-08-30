using DrWatson
@quickactivate 
using Revise
using Distributions
using LinearAlgebra
using Random
using SPM
using StatsBase



include(scriptsdir("cfg.jl"))

function simulate_control_chart(CH; target=mean, title="sims", statname="", h_up, ncond, f_tol, x_tol, maxrl, seed)
    # Precompilation
    bisectionCL(CH, h_up, nsims = 1, maxrl = 1)
    approximateBisectionCL(CH, nsims = 1, maxrl = 1)
    saCL(CH, maxiter = 1)
    time()
    # end precompilation

    if length(statname) == 0
        statname = string(typeof(STAT).name.wrapper)
    end

    Random.seed!(seed)
    nsims_vec = [1000, 10000]
    B_vec = [[100, 250, 500, 1000, 5000], [1000, 2500, 5000, 10000, 50000]]
    for i in 1:length(nsims_vec)
        nsims = nsims_vec[i] 
        fname = datadir(title, statname, "bisection", "sim-$(seed)-$(nsims).jld2")
        if !isfile(fname)
            t_bisec = time()
            bisectionCL!(CH, h_up, nsims = nsims, maxrl = maxrl, x_tol=x_tol, f_tol=f_tol)
            dt_bisec = time() - t_bisec
            RLs_bisec = zeros(ncond);
            for i in 1:ncond
                RLs_bisec[i] = run_sim(CH)
            end
            safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => target(RLs_bisec), "nsims" => nsims, "B" => nsims, "time" => dt_bisec))
        end

        RLs_approx = zeros(ncond);
        for B in B_vec[i]
            fname = datadir(title, statname, "approximate", "sim-$(seed)-$(nsims)-$(B).jld2")
            if !isfile(fname)
                t_approx = time()
                approximateBisectionCL!(CH, nsims=nsims, B = B, maxrl = maxrl, x_tol=x_tol, f_tol=f_tol)
                dt_approx = time() - t_approx
                for i in 1:ncond
                    RLs_approx[i] = run_sim(CH)
                end
                safesave(fname, Dict("h" => get_h(get_limit(CH)), "target" => target(RLs_approx), "nsims" => nsims, "B" => B, "time" => dt_approx))
            end
        end
    end
end


#############################################################
#                          LLCUSUM                          #
#############################################################

index_sim = parsed_args["index"]
seed = seeds[1] + index_sim
h_up = 100.0
ncond = 10_000
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

simulate_control_chart(CH, target=mean, h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)



#----- LLCUSUM median
NOM = QRL(200, 0.5)
STAT = LLCUSUM(1.0, x)
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2(MultinomialBootstrap(STAT), x)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_control_chart(CH, target=median, statname = "LLCUSUM-median", h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)


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

simulate_control_chart(CH, target=mean, h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)



#----- MEWMA median
NOM = QRL(200, 0.5)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_control_chart(CH, target=median, statname="MEWMA-median", h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)



#############################################################
#                       AMCUSUM 5 dim                       #
#############################################################
p = 5
seed = seeds[3] + index_sim
STAT = MCUSUM(0.25, p, 0.0, zeros(p))
LIM = OneSidedFixedLimit(5.0, true)
BOOT = Phase2Distribution(MvNormal(zeros(p), ones(p)))
NOM = ARL(200)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_control_chart(CH, target=mean, statname="AMCUSUM", h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)



#----- AMCUSUM median
NOM = QRL(200, 0.5)
CH = ControlChart(STAT, LIM, NOM, BOOT)

simulate_control_chart(CH, target=median, statname="AMCUSUM-median", h_up=h_up, ncond=ncond, f_tol=f_tol, x_tol=x_tol, maxrl=maxrl, seed=seed)
