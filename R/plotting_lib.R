library(ggplot2)
library(latex2exp)


#' @export
PlotInfluence <- function(influence_dfs, alpha_colname, alpha_max,
                          target_change=NULL,
                          include_y_zero=TRUE, sig_num_ses=qnorm(0.975)) {
  alpha_col <- sym(alpha_colname)
  base_vals <- SafeGetBaseVals(influence_dfs)

  influence_df <-
    bind_rows(influence_dfs$pos %>% mutate(sort="pos"),
              influence_dfs$neg %>% mutate(sort="neg")) %>%
    filter(!!alpha_col <= alpha_max)

  plot <- ggplot(influence_df, aes(x=!!alpha_col))
  if (!is.null(target_change)) {
    for (row in 1:nrow(target_change)) {
      this_alpha <- target_change[row, alpha_colname]
      if (!is.na(this_alpha) && (this_alpha < alpha_max)) {
        this_change <- target_change[row, "change"]
        plot <- plot + geom_vline(aes(xintercept=!!this_alpha,
                                      linetype=!!this_change))
      }
    }
  }
  if (include_y_zero) {
    plot <-
      plot +
      geom_line(aes(y=0.0), col="gray50")
  }
  plot <-
    plot +
    geom_line(aes(y=base_vals["beta"]), col="blue", lwd=2) +
    geom_ribbon(aes(
      ymin=beta_est - sig_num_ses * se_est,
      ymax=beta_est + sig_num_ses * se_est,
      group=sort),
      fill="blue", color=NA, alpha=0.1)

  xlab_name <-
    case_when(alpha_colname == "prop_removed" ~ "Proportion of points removed",
              alpha_colname == "num_removed" ~ "Number of points removed",
              TRUE ~ alpha_colname)
  plot <-
    plot +
    geom_line(aes(y=beta_est,
                  group=sort,
                  color="Observed"), lwd=2) +
    xlab(xlab_name) + ylab(TeX("$\\beta$")) +
    guides(color=FALSE,
           linetype=guide_legend(title="Change type"))

  return(plot)
}
