###################################################################
#
# These simple examples illustrate the use of zaminfluence.
# https://github.com/rgiordan/zaminfluence
# See the README.md file for installation instructions.

library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sandwich)
library(zaminfluence)
library(AER)

compare <- function(x, y) { return(max(abs(x - y))) }
check_equivalent  <- function(x, y) { stopifnot(compare(x, y) < 1e-8) }

n_obs <- 10000

set.seed(42)

git_repo_dir <- "/home/rgiordan/Documents/git_repos/zaminfluence"

#source(file.path(git_repo_dir, "zaminfluence/R/influence_lib.R"))


#############################
# Oridinary regression.

x_dim <- 3
beta_true <- 0.1 * runif(x_dim)
df <- GenerateRegressionData(n_obs, beta_true, num_groups=NULL)

# Fit a regression model.
x_names <- sprintf("x%d", 1:x_dim)
reg_form <- formula(sprintf("y ~ %s - 1", paste(x_names, collapse=" + ")))
reg_fit <- lm(data = df, formula = reg_form, x=TRUE, y=TRUE)

# Get influence.
reg_infl <- ComputeModelInfluence(reg_fit)
reg_infl <- AppendTargetRegressorInfluence(reg_infl, "x1")
reg_signals <-
    GetRegressionSignals(reg_infl$targets[["x1"]]) %>%
    RerunForTargetChanges(reg_infl)

# Make the reruns.
rerun_df <- 
    lapply(c("sign", "sig", "both"), 
           function(x) { reg_signals[[x]]$rerun_df }) %>%
    do.call(bind_rows, .)


parameter_infl <- reg_infl$targets[["x1"]]
names(parameter_infl)
sorting_qoi_name <- "beta"
qoi_infl <- parameter_infl[[sorting_qoi_name]]
qoi_infl %>% names()
qoi_infl$neg %>% names()
plot_num_dropped <- FALSE

##################

drop_inds <- parameter_infl[[qoi_for_sorting]]$neg$infl_inds[1:100]

PredictChange <- function(qoi_infl, drop_inds) {
    return(sum(qoi_infl$infl[drop_inds]))
}



GetSortedInfluenceDf <- function(parameter_infl, sorting_qoi_name, max_num_obs=Inf) {
    qoi_for_sorting <- parameter_infl[[sorting_qoi_name]]

    GetQOIDf <- function(infl_sign) {
        ordered_inds <- qoi_for_sorting[[infl_sign]]$infl_inds
        if (max_num_obs < length(ordered_inds)) {
            ordered_inds <- ordered_inds[1:max_num_obs]
        }
        qoi_df <- data.frame(
            num_dropped=c(0, 1:length(ordered_inds))) %>%
            mutate(prop_dropped=num_dropped / qoi_for_sorting[[infl_sign]]$num_obs,
                   sign=infl_sign)
        for (qoi_name in c("beta", "beta_mzse", "beta_pzse")) {
            base_value <- parameter_infl[[qoi_name]]$base_value
            infl_sorted <- parameter_infl[[qoi_name]]$infl[ordered_inds]
            qoi_df[[qoi_name]] <- base_value + cumsum(c(0, infl_sorted))
        }
        return(qoi_df)
    }
    
    qoi_df <-
        bind_rows(GetQOIDf("pos"), GetQOIDf("neg")) %>%
        mutate(sorted_by=sorting_qoi_name)
    
    return(qoi_df)
}


influence_df <- GetSortedInfluenceDf(parameter_infl, "beta")
head(influence_df)

alpha_colname <- "num_dropped"
alpha_max <- 300

alpha_col <- sym(alpha_colname)



plot_num_dropped <- FALSE
include_y_zero <- TRUE
sig_num_ses <- 2
apip_max <- 300



plot <- ggplot(influence_df %>% filter(alpha <= alpha_max), aes(x=alpha)) +
  
if (include_y_zero) {
    plot <-
        plot +
        geom_line(aes(y=0.0), col="gray50")
}

reg_signal <- reg_signals[["sign"]]  
plot <- PlotRegSignal(plot, 
                      reg_signals[["sign"]],
                      plot_num_dropped=plot_num_dropped,
                      apip_max=apip_max)
print(plot)



PlotInfluence <- function(influence_df,
                          plot_num_dropped=FALSE,
                          apip_max=NULL,
                          reg_signals=NULL,
                          include_y_zero=TRUE) {

    PlotRegSignal <- function(plot, reg_signal) {
        alpha_type <- if (plot_num_dropped) "n" else "prop"
        alpha <- reg_signal$apip[[alpha_type]]
        if (is.null(apip_max) || (!is.null(apip_max) && alpha <= apip_max)) {
            plot <- plot + geom_vline(aes(xintercept=!!alpha,
                                          linetype=!!reg_signal$change))
            if (!is.null(reg_signal$rerun_df)) {
                rerun_df <- reg_signal$rerun_df
                stopifnot(nrow(rerun_df) == 1)
                plot <- 
                    plot + 
                    geom_errorbar(aes(
                        x=!!alpha, 
                        ymin=beta_mzse_refit,
                        ymax=beta_pzse_refit),
                        data=rerun_df,
                        width=errorbar_width,
                        lwd=1.5) +
                    geom_point(aes(x=!!alpha, y=betahat_refit),
                               data=rerun_df,
                               shape=8)
            }
        }
        return(plot)
    }

    influence_df$alpha <- 
        if (plot_num_dropped) influence_df$num_dropped else 
            influence_df$prop_dropped
    influence_df <- filter(influence_df, alpha  <= apip_max)
    errorbar_width <- diff(range(influence_df$alpha)) / 50
    
    plot <- ggplot(influence_df, aes(x=alpha)) 
    if (include_y_zero) {
        plot <-
            plot +
            geom_line(aes(y=0.0), col="gray50")
    }
    if (!is.null(reg_signals)) {
        for (reg_signal in reg_signals) {
            plot <- PlotRegSignal(plot,  reg_signal)
        }
        plot <- plot + guides(linetype=guide_legend(title="Change type"))
    }
    
    base_beta <- filter(influence_df, alpha == 0) %>% pull("beta") %>% unique()
    stopifnot(length(base_beta) == 1)
    plot <-
        plot +
        geom_line(aes(y=!!base_beta), col="blue", lwd=2) +
        geom_ribbon(aes(
            ymin=beta_mzse,
            ymax=beta_pzse,
            group=sign),
            fill="blue", color=NA, alpha=0.1) +
        geom_line(aes(y=beta, group=sign), lwd=2)
    
    xlab_name <- if (plot_num_dropped) 
        "Number of points removed" else "Proportion of points removed"
    plot <- plot + guides(color="none") + xlab(xlab_name)
    
    return(plot)
}

#debug(PlotInfluence)
PlotInfluence(influence_df, reg_signals=list(reg_signals[["sign"]]), apip_max=0.05)

