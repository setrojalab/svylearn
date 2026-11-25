.need <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the optional package '", package, "' to use this function.", call. = FALSE)
}
.positive <- function(w, n, zero = FALSE) {
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w)) ||
      any(if (zero) w < 0 else w <= 0) || sum(w) <= 0)
    stop("Weights must be finite, correctly sized and ",
         if (zero) "nonnegative with positive sum." else "strictly positive.", call. = FALSE)
  w
}
.binary <- function(y) {
  if (!(is.numeric(y) || is.logical(y)) || !length(y) ||
      anyNA(y) || any(!y %in% c(0, 1)))
    stop("Outcome must be a nonmissing numeric or logical 0/1 vector.", call. = FALSE)
  as.numeric(y)
}
.prob <- function(p, n) {
  if (!is.numeric(p) || length(p) != n || any(!is.finite(p)) || any(p < 0 | p > 1))
    stop("Probabilities must be finite, in [0, 1], and match the outcome length.", call. = FALSE)
  p
}
.clip <- function(p) pmax(1e-15, pmin(1 - 1e-15, p))
.seed <- function(seed, code) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}
.cols <- function(x) if (inherits(x, "formula")) all.vars(x) else x
.wquantile <- function(x, w, probs) {
  ok <- is.finite(x) & w > 0
  if (!any(ok)) stop("A numeric predictor has no observed training values.", call. = FALSE)
  x <- x[ok]; w <- w[ok]; o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}
.key <- function(data, cols) {
  if (!length(cols)) return(rep("1", nrow(data)))
  # Length-prefixed components avoid collisions from delimiters in identifiers.
  do.call(paste0, lapply(data[cols], function(v) {
    v <- as.character(v); paste0(nchar(v), ":", v, ";")
  }))
}
