#' Learn weighted preprocessing using training observations only
#' @param data Training data frame of predictors only.
#' @param w Training weights.
#' @param scale Standardize numeric columns using weighted training moments.
#' @param winsor Optional pair of probability bounds for numeric winsorization.
#' @param spline_df Natural spline degrees of freedom, zero to disable.
#' @param pca_rank Optional maximum number of weighted principal components after
#'   encoding. Centering, singular vectors and retained rank are learned on training
#'   data only. Uses a thin SVD, without constructing a predictor covariance matrix.
#' @return A fitted svy_recipe object. Apply it with bake_svy_recipe.
#' @export
svy_recipe <- function(data, w = rep(1, nrow(data)), scale = TRUE,
                       winsor = NULL, spline_df = 0L, pca_rank = NULL) {
  if (!is.data.frame(data) || !ncol(data)) stop("Supply a predictor data frame.")
  .positive(w, nrow(data))
  if (!is.null(winsor) && (length(winsor) != 2L || anyNA(winsor) ||
      winsor[1] < 0 || winsor[2] > 1 || winsor[1] >= winsor[2])) stop("Invalid winsor bounds.")
  if (length(spline_df) != 1L || spline_df < 0L || spline_df != as.integer(spline_df)) stop("Invalid spline_df.")
  specs <- lapply(data, function(x) {
    if (is.numeric(x)) {
      med <- .wquantile(x, w, .5); x[!is.finite(x)] <- med
      lim <- if (is.null(winsor)) c(-Inf, Inf) else .wquantile(x, w, winsor)
      x <- pmin(lim[2], pmax(lim[1], x))
      mu <- if (scale) sum(w * x) / sum(w) else 0
      sd <- if (scale) sqrt(sum(w * (x - mu)^2) / sum(w)) else 1
      if (sd < sqrt(.Machine$double.eps)) sd <- 1
      z <- (x - mu) / sd
      basis <- if (spline_df > 0L && length(unique(z)) > spline_df)
        splines::ns(z, df = spline_df) else NULL
      list(type = "numeric", median = med, bounds = lim, mean = mu, sd = sd, basis = basis)
    } else {
      x <- as.character(x); ok <- !is.na(x)
      if (!any(ok)) stop("A categorical predictor has no observed training values.")
      freq <- tapply(w[ok], x[ok], sum)
      mode <- names(freq)[which.max(freq)]; lev <- sort(unique(x[ok]))
      # Novel categories map to the training mode, with no validation-derived levels.
      list(type = "factor", mode = mode, levels = c(mode, setdiff(lev, mode)))
    }
  })
  rec <- structure(list(specs = specs, names = names(data)), class = "svy_recipe")
  if (!is.null(pca_rank)) {
    if (length(pca_rank) != 1L || !is.finite(pca_rank) || pca_rank < 1L || pca_rank != as.integer(pca_rank))
      stop("pca_rank must be a positive integer.")
    x <- bake_svy_recipe(rec, data); ww <- w / sum(w)
    center <- colSums(x * ww); z <- sweep(x, 2L, center)
    max_rank <- min(pca_rank, ncol(z), nrow(z) - 1L)
    if (max_rank < 1L) stop("PCA requires at least two training rows.")
    s <- base::svd(z * sqrt(ww), nu = 0L, nv = max_rank)
    rank <- min(max_rank, sum(s$d > max(s$d) * sqrt(.Machine$double.eps)))
    if (rank < 1L) stop("PCA has no nonconstant training directions.")
    rotation <- s$v[, seq_len(rank), drop = FALSE]
    sign <- vapply(seq_len(rank), function(j) sign(rotation[which.max(abs(rotation[,j])),j]), numeric(1))
    rotation <- sweep(rotation, 2L, sign, `*`)
    rownames(rotation) <- colnames(x); colnames(rotation) <- paste0("PC", seq_len(rank))
    rec$pca <- list(center = center, rotation = rotation, rank = rank,
      variance = s$d[seq_len(rank)]^2, explained_fraction = s$d[seq_len(rank)]^2 / sum(s$d^2))
  }
  rec
}

#' Apply frozen training-fold preprocessing
#' @param recipe A fitted svy_recipe.
#' @param newdata Data frame with the required predictor columns.
#' @return Numeric design matrix without an intercept.
#' @export
bake_svy_recipe <- function(recipe, newdata) {
  if (!all(recipe$names %in% names(newdata))) stop("Required predictors missing from newdata.")
  ans <- lapply(recipe$names, function(nm) {
    s <- recipe$specs[[nm]]; x <- newdata[[nm]]
    if (s$type == "numeric") {
      if (!is.numeric(x)) stop("Numeric predictor changed type: ", nm)
      x[!is.finite(x)] <- s$median
      z <- (pmin(s$bounds[2], pmax(s$bounds[1], x)) - s$mean) / s$sd
      m <- if (is.null(s$basis)) matrix(z, ncol = 1L) else stats::predict(s$basis, z)
      colnames(m) <- if (ncol(m) == 1L) nm else paste0(nm, "_s", seq_len(ncol(m)))
      m
    } else {
      x <- as.character(x); x[is.na(x) | !x %in% s$levels] <- s$mode
      lev <- s$levels[-1L]
      m <- matrix(0, nrow(newdata), length(lev))
      for (j in seq_along(lev)) m[, j] <- as.numeric(x == lev[j])
      if (length(lev)) colnames(m) <- paste0(nm, "=", lev)
      m
    }
  })
  m <- do.call(cbind, ans)
  if (!ncol(m)) stop("No nonconstant encoded predictors.")
  if (!is.null(recipe$pca)) m <- sweep(m, 2L, recipe$pca$center) %*% recipe$pca$rotation
  m
}
