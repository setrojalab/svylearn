#' Export prediction metrics to a text table
#' @param x A data frame or svy_learn_fit.
#' @param file Output file path.
#' @param format csv, markdown or latex.
#' @param digits Printed decimal places.
#' @return The output path invisibly.
#' @export
svy_export <- function(x, file, format = "csv", digits = 4L) {
  format <- match.arg(format, c("csv", "markdown", "latex"))
  if (inherits(x, "svy_learn_fit")) x <- data.frame(metric = colnames(x$metrics),
    mean = colMeans(x$metrics, na.rm = TRUE), sd = apply(x$metrics, 2, stats::sd), row.names = NULL)
  if (!is.data.frame(x)) stop("x must be a data frame or fitted object.")
  if (format == "csv") {utils::write.csv(x, file, row.names = FALSE); return(invisible(file))}
  x[] <- lapply(x, function(z) if (is.numeric(z)) formatC(z, digits = digits, format = "f") else as.character(z))
  if (format == "markdown") {
    escape <- function(z) gsub("|", "\\|", z, fixed = TRUE)
    x[] <- lapply(x, escape)
    rows <- apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
    lines <- c(paste0("| ", paste(escape(names(x)), collapse = " | "), " |"),
      paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"), rows)
  } else {
    escape <- function(z) {
      z <- gsub("\\", "\\textbackslash{}", z, fixed = TRUE)
      for (s in c("_", "%", "&", "#", "$")) z <- gsub(s, paste0("\\", s), z, fixed = TRUE)
      z
    }
    x[] <- lapply(x, escape)
    lines <- c(paste0("\\begin{tabular}{", paste(rep("l", ncol(x)), collapse = ""), "}"),
      paste0(paste(escape(names(x)), collapse = " & "), " \\\\"), "\\hline",
      apply(x, 1, function(z) paste0(paste(z, collapse = " & "), " \\\\")), "\\end{tabular}")
  }
  writeLines(lines, file, useBytes = TRUE); invisible(file)
}
