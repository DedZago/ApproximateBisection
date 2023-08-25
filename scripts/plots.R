require(pacman)
pacman::p_load("ggplot2", "mvbutils", "RColorBrewer")
pacman::p_load("knitr", "kableExtra", "flextable", "officer", "tidyverse")
theme_set(theme_bw())
setwd("/home/dzago/ApproximateBisection/data/output/sims")

LOADDIR = getwd()
for(statistic in list.dirs(LOADDIR, recursive = FALSE)){
    setwd(statistic)
    cat("-----", statistic, "-----", "\n")
    df = read.csv("results_rounded.csv")
    print(df)
    cat("\n\n")
    setwd("..")
}