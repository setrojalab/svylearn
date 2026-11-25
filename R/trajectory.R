#' Construct time-eligible subject trajectory features
#' @param data Long-format data frame.
#' @param subject,wave Subject and finite numeric time column names.
#' @param variables Numeric measurement columns.
#' @param cutoff Last eligible numeric time, specified before evaluating outcomes.
#' @return One row per subject with observed count, mean, last observation and
#'   least-squares time slope for each measurement. Subjects with no eligible
#'   observations have missing features; impute using training-fold recipes.
#'   No survey weights are aggregated or inferred. Future measurements never enter
#'   the summaries. Slopes require at least two distinct observed times.
#' @export
svy_trajectory <- function(data, subject, wave, variables, cutoff) {
  if (!is.data.frame(data) || !nrow(data) || length(subject) != 1L || length(wave) != 1L ||
      !length(variables) || anyDuplicated(variables) ||
      !all(c(subject, wave, variables) %in% names(data)) || any(variables %in% c(subject,wave)))
    stop("Invalid trajectory columns.")
  if (anyNA(data[[subject]]) || !is.numeric(data[[wave]]) || any(!is.finite(data[[wave]])) ||
      anyDuplicated(data[c(subject,wave)])) stop("Finite times, observed subjects and unique subject-wave rows are required.")
  if (!all(vapply(data[variables],is.numeric,logical(1))) || length(cutoff) != 1L || !is.finite(cutoff))
    stop("Numeric measurements and a finite prespecified cutoff are required.")
  subjects <- unique(data[[subject]])
  rows <- split(seq_len(nrow(data)),match(data[[subject]],subjects))
  rows <- rows[as.character(seq_along(subjects))]
  out <- data.frame(subjects,check.names=FALSE); names(out) <- subject
  for (variable in variables) {
    summaries <- t(vapply(rows,function(i) {
      i <- i[data[[wave]][i] <= cutoff & is.finite(data[[variable]][i])]
      if (!length(i)) return(c(n=0,mean=NA_real_,last=NA_real_,slope=NA_real_))
      time <- data[[wave]][i]; value <- data[[variable]][i]
      slope <- if (length(i)>1L) sum((time-mean(time))*(value-mean(value)))/sum((time-mean(time))^2) else NA_real_
      c(n=length(i),mean=mean(value),last=value[which.max(time)],slope=slope)
    },numeric(4)))
    colnames(summaries) <- paste(variable,colnames(summaries),sep="_")
    out <- cbind(out,summaries)
  }
  attr(out,"cutoff") <- cutoff; out
}
