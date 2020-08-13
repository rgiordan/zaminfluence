library(dplyr)


PairInfluenceScoresForSign <- function(influence_df, assignment_col, sort_col,
                                       influence_cols, change_sign,
                                       group_cols=NULL,
                                       level0=0, level1=1) {
  change_sign <- sign(change_sign)
  if (change_sign == 0) {
    stop("change_sign must be nonzero.")
  }
  if (!(sort_col %in% names(influence_df))) {
    stop(sprintf("sort_col %s is not a column in influence_df", sort_col))
  }
  if (!(assignment_col %in% names(influence_df))) {
    stop(sprintf("assignment_col %s is not a column in influence_df",
                 assignment_col))
  }
  for (influence_col in influence_cols) {
    influence_col <- paste0(influence_col, "_grad")
    if (!(influence_col %in% names(influence_df))) {
      stop(sprintf("influence_col %s is not a column in influence_df",
                   influence_col))
    }
  }
  if (!is.null(group_cols)) {
    for (group_col in group_cols) {
      if (!(group_col %in% names(influence_df))) {
        stop(sprintf("group_col %s is not a column in influence_df", group_col))
      }
    }
  }

  n_obs <- attr(influence_df, "n_obs")

  level0_rows <- influence_df[[assignment_col]] == level0
  level1_rows <- influence_df[[assignment_col]] == level1
  if (sum(level0_rows) == 0) {
    stop("Level 0 has no rows.")
  }
  if (sum(level1_rows) == 0) {
    stop("Level 1 has no rows.")
  }
  active_rows <- level0_rows || level1_rows
  assignment_sym <- sym(assignment_col)
  sort_sym <- sym(sort_col)
  obs_per_row_col <- attr(influence_df, "obs_per_row_col")

  ranked_df <-
    influence_df[active_rows, ] %>%
    group_by(!!assignment_sym)

  # Optionally add other grouping columns and rank within them....
  if (!is.null(group_cols)) {
    for (group_col in group_cols) {
      group_sym <- sym(group_col)
      ranked_df <- group_by(ranked_df, !!group_sym, .add=TRUE)
    }
  }

  ranked_df <-
    ranked_df %>%
    mutate(rank=as.integer(rank(change_sign * !!sort_sym,
                                ties.method="random")))

  final_paired_df <- inner_join(
      filter(ranked_df, !!assignment_sym == !!level0),
      filter(ranked_df, !!assignment_sym == !!level1),
      by=c("rank", group_cols),
      suffix=c("_0", "_1"))

  final_paired_df[[sort_col]] <-
    final_paired_df[[paste0(sort_col, "_1")]] +
    final_paired_df[[paste0(sort_col, "_0")]]

  for (influence_col in influence_cols) {
    influence_col <- paste0(influence_col, "_grad")
    final_paired_df[[influence_col]] <-
      final_paired_df[[paste0(influence_col, "_1")]] +
      final_paired_df[[paste0(influence_col, "_0")]]
  }

  # Count the number removed.
  final_paired_df[[obs_per_row_col]] <-
    final_paired_df[[paste0(obs_per_row_col, "_1")]] +
    final_paired_df[[paste0(obs_per_row_col, "_0")]]

  final_paired_df <-
    final_paired_df %>%
    ungroup() %>%   # Necessary to get aggregation correct
    arrange(change_sign * !!sort_sym) %>%
    mutate(num_removed=cumsum(get(obs_per_row_col))) %>%
    mutate(prop_removed=num_removed / n_obs) %>%
    PrependZeroRemoved(c(sort_col, paste0(influence_cols, "_grad")))

  base_vals <- attr(influence_df, "base_vals")
  final_paired_df <-
    final_paired_df %>%
    AccumulateInfluenceDf(
      influence_cols=paste0(influence_cols, "_grad"),
      base_vals=base_vals[influence_cols],
      change_colnames=paste0(influence_cols, "_change"),
      estimate_colnames=paste0(influence_cols, "_est")) %>%
      CopyGradAttributes(influence_df)

  attr(final_paired_df, "sort_col") <- sort_col
  attr(final_paired_df, "change_sign") <- change_sign
  attr(final_paired_df, "influence_cols") <- influence_cols
  attr(final_paired_df, "group_cols") <- group_cols
  attr(final_paired_df, "assignment_col") <- assignment_col
  attr(final_paired_df, "level0") <- level0
  attr(final_paired_df, "level1") <- level1
  attr(final_paired_df, "obs_per_row_col") <- obs_per_row_col
  attr(final_paired_df, "data_row_cols") <-
    c(paste0(attr(influence_df, "data_row_cols"), "_0"),
      paste0(attr(influence_df, "data_row_cols"), "_1"))

  return(final_paired_df)
}


#' Process influence dataframes to require observations to be removed in pairs.
#' @param influence_df The output of GetTargetRegressorGrads with assignment
#' information
#' @param assignment_col A string containing the column name of the assignment
#' @param influence_cols Optional custom influence column names
#' @param group_cols Optional.  If you wish to make sure that additional
#' attributes match in removed pairs, you can specify the columns to match on
#' by passing the column names in as group_cols.
#' @param level0 The first value of assignment.  Each pair will have one
#' with this assignment level.
#' @param level1 The second value of assignment.  Each pair will also have one
#' with this assignment level.
#' @return A list of influence dataframes analogous to that returned by
#' SortAndAccumulate, but with paired, grouped observations.
#' @export
PairInfluenceScores <- function(influence_df,
                                assignment_col,
                                influence_cols=NULL,
                                group_cols=NULL,
                                level0=0, level1=1) {
  if (is.null(influence_cols)) {
      influence_cols <- names(attr(influence_df, "base_vals"))
  }
  LocalPairInfluenceScoresForSign <- function(sort_col, change_sign) {
    PairInfluenceScoresForSign(
      influence_df=influence_df,
      assignment_col=assignment_col,
      influence_cols=influence_cols,
      sort_col=sort_col,
      group_cols=group_cols,
      change_sign=change_sign)
  }
  result <- list(
    sign=list(
      pos=LocalPairInfluenceScoresForSign("beta_grad", change_sign=1),
      neg=LocalPairInfluenceScoresForSign("beta_grad", change_sign=-1)),
    sig=list(
      pos=LocalPairInfluenceScoresForSign("beta_pzse_grad", change_sign=1),
      neg=LocalPairInfluenceScoresForSign("beta_mzse_grad", change_sign=-1))
  )
  return(result)
}
