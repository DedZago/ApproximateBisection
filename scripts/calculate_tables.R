require(pacman)
pacman::p_load("ggplot2", "mvbutils", "RColorBrewer", "data.table")
pacman::p_load("knitr", "kableExtra", "flextable", "officer", "tidyverse")
theme_set(theme_bw())

################### Comparison with standard bisection ###################
setwd("/home/dzago/ApproximateBisection/data/output/sims")
LOADDIR = getwd()
df_final = data.frame()
statistics = list.dirs(LOADDIR, recursive = FALSE)

for(i in 1:lens){
    setwd(statistics[i])
    cat("-----", statistics[i], "-----", "\n")
    df = read.csv("results_rounded.csv")
    rn = df[,1]
    rn = str_replace(rn, "bisection", "Bisection")
    rn = str_replace(rn, "approximate", "BA-Bisection")
    df[,1] = rn
    df = df[, c(1,2,3,4,7,5,8,6,9)]
    cn = colnames(df)
    print(df)
    linesep = c(rep("", 5), "\\addlinespace", rep("",5))
    cat("\n\n")
    tab <- kable(df, format="latex", booktabs=TRUE, digits = 3, row.names=FALSE, escape=FALSE, align='c',
            linesep = linesep, caption = "", label="results") %>%
            add_header_above(c("", "", "", "time (s)" = 2, "target" = 2, "h" = 2)) %>%
            kable_styling(latex_options = "hold_position")

    writeLines(tab, paste0(statistics[i], "/results.tex"))
    setwd("..")
}

################### Comparison with SACL ###################
setwd("/home/dzago/ApproximateBisection/data/output/sims-sacl")
LOADDIR = getwd()
df_final = data.frame()
statistics = list.dirs(LOADDIR, recursive = FALSE)
statistics = statistics[c(5,6,1,2,3,4)]
lens = length(statistics)

df_list = vector(mode="list", length = lens)
for(i in 1:lens){
    setwd(statistics[i])
    cat("-----", statistics[i], "-----", "\n")
    df = read.csv("results_rounded.csv")
    rn = df[,1]
    cn = colnames(df)
    df_t = t(df)
    rownames(df_t) = cn
    colnames(df_t) = rn
    df_t = df_t[-1,]
    N = NROW(df_t)
    df_tab = cbind(df_t[1:(N/2),2], df_t[(N/2+1):N,2])
    df_tab = cbind(df_tab, df_t[1:(N/2),1], df_t[(N/2+1):N,1])
    df_tab = as.data.frame(df_tab); df_tab = apply(df_tab, 2, as.numeric)
    cn_tab = cn[2:(N/2+1)]
    cn_tab = str_replace(cn_tab, "_median", "")
    cn_tab = str_replace(cn_tab, "_median", "")
    rownames(df_tab) = cn_tab
    colnames(df_tab) = rep(c("median", "iqr"), 2)
    print(df_tab)
    cat("\n\n")
    setwd("..")
    df_list[[i]] = df_tab
}

# df_all = do.call("rbind", df_list)
row_counts = sapply(df_list, NROW)
stat_counts = rep(NA, 3)
# stat_counts[1] = row_counts[1] + row_counts[2]
# stat_counts[2] = row_counts[3] + row_counts[4]
# stat_counts[3] = row_counts[5] + row_counts[6]
names(row_counts) = rep(c("ARL", "MRL"), length(row_counts)/2)
df_all = rbind(cbind(df_list[[1]], df_list[[2]]), cbind(df_list[[3]], df_list[[4]]), cbind(df_list[[5]], df_list[[6]]))
stat_counts[1] = row_counts[1]
stat_counts[2] = row_counts[3]
stat_counts[3] = row_counts[5]
names(stat_counts) = c("Multiple EWMA", "Multiple LLCUSUM", "T2 and multiple MCUSUM")

linesep = do.call("c", sapply(stat_counts, function(n) c(rep('', n-1), "\\addlinespace")))
tab <- kable(df_all, format="latex", booktabs=TRUE, digits = 3, row.names=TRUE, escape=FALSE, align='c',
        linesep = linesep, caption = "", label="SA-results") %>%
        add_header_above(c("", "BA-Bisection" = 2, "SA"=2, "BA-Bisection" = 2, "SA" = 2)) %>%
        add_header_above(c("", "ARL" = 4, "MRL" = 4)) %>%
        pack_rows(index = stat_counts) %>%
        kable_styling(latex_options = "hold_position")


################### RSADA 4 CORES ###################
setwd("/home/dzago/ApproximateBisection/data/output/sims-parallel")
LOADDIR = getwd()
df_final = data.frame()
statistics = list.dirs(LOADDIR, recursive = FALSE)
lens = length(statistics)

df_list = vector(mode="list", length = lens)
for(i in 1:lens){
    setwd(statistics[i])
    cat("-----", statistics[i], "-----", "\n")
    df = read.csv("results_rounded.csv")
    rn = df[,1]
    cn = colnames(df)
    df_t = t(df)
    rownames(df_t) = cn
    colnames(df_t) = rn
    df_t = df_t[-1,,drop=FALSE]
    N = NROW(df_t)
    df_tab = cbind(df_t[1:(N/2),1], df_t[(N/2+1):N,1])
    df_tab = as.data.frame(df_tab); df_tab = apply(df_tab, 2, as.numeric)
    cn_tab = cn[2:(N/2+1)]
    cn_tab = str_replace(cn_tab, "_median", "")
    cn_tab = str_replace(cn_tab, "_median", "")
    rownames(df_tab) = cn_tab
    colnames(df_tab) = rep(c("median", "iqr"), 1)
    print(df_tab)
    cat("\n\n")
    setwd("..")
    df_list[[i]] = df_tab
}

# df_all = do.call("rbind", df_list)
row_counts = sapply(df_list, NROW)
stat_counts = rep(NA, 3)
# stat_counts[1] = row_counts[1] + row_counts[2]
# stat_counts[2] = row_counts[3] + row_counts[4]
# stat_counts[3] = row_counts[5] + row_counts[6]
names(row_counts) = rep(c("ARL", "MRL"), length(row_counts)/2)
df_all = rbind(cbind(df_list[[1]], df_list[[2]]))

tab <- kable(df_all, format="latex", booktabs=TRUE, digits = 3, row.names=TRUE, escape=FALSE, align='c',
        linesep = "", caption = "", label="RSADA-results") %>%
        add_header_above(c("", "ARL" = 2, "MRL" = 2)) %>%
        kable_styling(latex_options = "hold_position")


################### RSADA 8 CORES ###################
setwd("/home/dzago/ApproximateBisection/data/output/sims-parallel-8cores")
LOADDIR = getwd()
df_final = data.frame()
statistics = list.dirs(LOADDIR, recursive = FALSE)
lens = length(statistics)

df_list = vector(mode="list", length = lens)
for(i in 1:lens){
    setwd(statistics[i])
    cat("-----", statistics[i], "-----", "\n")
    df = read.csv("results_rounded.csv")
    rn = df[,1]
    cn = colnames(df)
    df_t = t(df)
    rownames(df_t) = cn
    colnames(df_t) = rn
    df_t = df_t[-1,,drop=FALSE]
    N = NROW(df_t)
    df_tab = cbind(df_t[1:(N/2),1], df_t[(N/2+1):N,1])
    df_tab = as.data.frame(df_tab); df_tab = apply(df_tab, 2, as.numeric)
    cn_tab = cn[2:(N/2+1)]
    cn_tab = str_replace(cn_tab, "_median", "")
    cn_tab = str_replace(cn_tab, "_median", "")
    rownames(df_tab) = cn_tab
    colnames(df_tab) = rep(c("median", "iqr"), 1)
    print(df_tab)
    cat("\n\n")
    setwd("..")
    df_list[[i]] = df_tab
}

row_counts = sapply(df_list, NROW)
stat_counts = rep(NA, 3)
names(row_counts) = rep(c("ARL", "MRL"), length(row_counts)/2)
df_all = rbind(cbind(df_list[[1]], df_list[[2]]))

tab <- kable(df_all, format="latex", booktabs=TRUE, digits = 3, row.names=TRUE, escape=FALSE, align='c',
        linesep = "", caption = "", label="RSADA-results-8cores") %>%
        add_header_above(c("", "ARL" = 2, "MRL" = 2)) %>%
        kable_styling(latex_options = "hold_position")

################### RSADA SEQUENTIAL ###################
setwd("/home/dzago/ApproximateBisection/data/output/sims-nonparallel")
LOADDIR = getwd()
df_final = data.frame()
statistics = list.dirs(LOADDIR, recursive = FALSE)
lens = length(statistics)

df_list = vector(mode="list", length = lens)
for(i in 1:lens){
    setwd(statistics[i])
    cat("-----", statistics[i], "-----", "\n")
    df = read.csv("results_rounded.csv")
    rn = df[,1]
    cn = colnames(df)
    df_t = t(df)
    rownames(df_t) = cn
    colnames(df_t) = rn
    df_t = df_t[-1,,drop=FALSE]
    N = NROW(df_t)
    df_tab = cbind(df_t[1:(N/2),1], df_t[(N/2+1):N,1])
    df_tab = as.data.frame(df_tab); df_tab = apply(df_tab, 2, as.numeric)
    cn_tab = cn[2:(N/2+1)]
    cn_tab = str_replace(cn_tab, "_median", "")
    cn_tab = str_replace(cn_tab, "_median", "")
    rownames(df_tab) = cn_tab
    colnames(df_tab) = rep(c("median", "iqr"), 1)
    print(df_tab)
    cat("\n\n")
    setwd("..")
    df_list[[i]] = df_tab
}

row_counts = sapply(df_list, NROW)
stat_counts = rep(NA, 3)
names(row_counts) = rep(c("ARL", "MRL"), length(row_counts)/2)
df_all = rbind(cbind(df_list[[1]], df_list[[2]]))

tab <- kable(df_all, format="latex", booktabs=TRUE, digits = 3, row.names=TRUE, escape=FALSE, align='c',
        linesep = "", caption = "", label="RSADA-results-nonparallel") %>%
        add_header_above(c("", "ARL" = 2, "MRL" = 2)) %>%
        kable_styling(latex_options = "hold_position")


