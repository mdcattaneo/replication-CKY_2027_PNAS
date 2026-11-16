# Streaming Monte Carlo moments. Tree fitting and data generation stay in the
# experiment entry points. Array order is grid point, depth, then method.
cky_mc_init <- function(grid, depths, methods) {
  ng <- nrow(grid)
  stopifnot(ng > 0L, length(depths) > 0L, length(methods) > 0L,
            !anyDuplicated(depths), !anyDuplicated(methods))
  pairs <- list()
  add_pair <- function(comparison, method_a, depth_a, method_b, depth_b) {
    pairs[[length(pairs) + 1L]] <<- data.frame(
      comparison = comparison, method_a = method_a, depth_a = depth_a,
      method_b = method_b, depth_b = depth_b, stringsAsFactors = FALSE)
  }
  if (length(depths) > 1L) {
    for (method in methods) for (k in 2:length(depths)) {
      add_pair("adjacent_depth", method, depths[k - 1L], method, depths[k])
    }
  }
  for (rule in c("IPW", "DIM", "SSE")) {
    a <- paste0("NSS-", rule)
    b <- paste0("HON-", rule)
    if (all(c(a, b) %in% methods)) for (depth in depths) {
      add_pair("HON_minus_NSS", a, depth, b, depth)
    }
  }
  pair_info <- if (length(pairs)) do.call(rbind, pairs) else data.frame(
    comparison = character(), method_a = character(), depth_a = integer(),
    method_b = character(), depth_b = integer())
  indices <- function(method, depth) {
    seq_len(ng) + ng * (match(depth, depths) - 1L +
                        length(depths) * (match(method, methods) - 1L))
  }
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  state$grid <- grid
  state$depths <- depths
  state$methods <- methods
  state$shape <- c(ng, length(depths), length(methods))
  state$error_mean <- state$error_m2 <- state$loss_mean <- state$loss_m2 <-
    numeric(prod(state$shape))
  state$pair_info <- pair_info
  state$pair_a <- unlist(Map(indices, pair_info$method_a, pair_info$depth_a), use.names = FALSE)
  state$pair_b <- unlist(Map(indices, pair_info$method_b, pair_info$depth_b), use.names = FALSE)
  state$pair_mean <- state$pair_m2 <- state$pair_comoment <- numeric(length(state$pair_a))
  state
}

cky_mc_update <- function(state, prediction, target) {
  if (!identical(as.integer(dim(prediction)), as.integer(state$shape)) ||
      length(target) != 1L || !is.finite(target) || any(!is.finite(prediction))) {
    stop("Invalid prediction shape, values or target for Monte Carlo accumulation")
  }
  error <- as.double(prediction) - target
  loss <- error * error
  if (any(!is.finite(loss))) stop("Nonfinite squared error")
  n <- state$n + 1L
  de <- error - state$error_mean
  state$error_mean <- state$error_mean + de / n
  state$error_m2 <- state$error_m2 + de * (error - state$error_mean)
  dl <- loss - state$loss_mean
  state$loss_mean <- state$loss_mean + dl / n
  state$loss_m2 <- state$loss_m2 + dl * (loss - state$loss_mean)
  a <- state$pair_a
  b <- state$pair_b
  difference <- loss[b] - loss[a]
  dd <- difference - state$pair_mean
  state$pair_mean <- state$pair_mean + dd / n
  state$pair_m2 <- state$pair_m2 + dd * (difference - state$pair_mean)
  state$pair_comoment <- state$pair_comoment + dl[a] * (loss[b] - state$loss_mean[b])
  state$n <- n
  invisible(state)
}

cky_mc_nonnegative <- function(value, scale = abs(value)) {
  # Allow roundoff in differences of highly correlated estimated moments.
  if (any(value < -1e-10 * pmax(1, scale), na.rm = TRUE)) {
    stop("Monte Carlo covariance produced a materially negative variance")
  }
  pmax(0, value)
}

cky_mc_pointwise <- function(state) {
  if (state$n < 1L) stop("No Monte Carlo replications accumulated")
  ng <- state$shape[1L]
  nk <- length(state$depths)
  nm <- length(state$methods)
  variance_error <- if (state$n > 1L) state$error_m2 / (state$n - 1L) else rep(NA_real_, length(state$error_m2))
  variance_loss <- if (state$n > 1L) state$loss_m2 / (state$n - 1L) else rep(NA_real_, length(state$loss_m2))
  mse <- state$loss_mean
  rmse <- sqrt(mse)
  mcse_mse <- sqrt(cky_mc_nonnegative(variance_loss) / state$n)
  mcse_rmse <- rep(NA_real_, length(rmse))
  positive <- rmse > 0
  mcse_rmse[positive] <- mcse_mse[positive] / (2 * rmse[positive])
  data.frame(
    state$grid[rep(seq_len(ng), nk * nm), , drop = FALSE],
    method = rep(state$methods, each = ng * nk),
    K = rep(rep(state$depths, each = ng), nm), replications = state$n,
    mean_error = state$error_mean,
    MCSE_mean_error = sqrt(cky_mc_nonnegative(variance_error) / state$n),
    error_variance = variance_error, MSE = mse, squared_loss_variance = variance_loss,
    RMSE = rmse, MCSE_MSE = mcse_mse, MCSE_RMSE = mcse_rmse,
    relative_MCSE_RMSE = mcse_rmse / rmse, row.names = NULL)
}

cky_mc_pairs <- function(state) {
  if (state$n < 1L) stop("No Monte Carlo replications accumulated")
  ng <- state$shape[1L]
  np <- nrow(state$pair_info)
  a <- state$pair_a
  b <- state$pair_b
  variance <- if (state$n > 1L) state$loss_m2 / (state$n - 1L) else rep(NA_real_, length(state$loss_m2))
  covariance <- if (state$n > 1L) state$pair_comoment / (state$n - 1L) else rep(NA_real_, length(a))
  variance_difference <- if (state$n > 1L) state$pair_m2 / (state$n - 1L) else rep(NA_real_, length(a))
  rmse_a <- sqrt(state$loss_mean[a])
  rmse_b <- sqrt(state$loss_mean[b])
  v_rmse <- rep(NA_real_, length(a))
  positive <- rmse_a > 0 & rmse_b > 0
  va <- variance[a[positive]] / (4 * state$loss_mean[a[positive]])
  vb <- variance[b[positive]] / (4 * state$loss_mean[b[positive]])
  cross <- covariance[positive] / (2 * rmse_a[positive] * rmse_b[positive])
  v_rmse[positive] <- cky_mc_nonnegative(va + vb - cross, abs(va) + abs(vb) + abs(cross))
  data.frame(
    state$grid[rep(seq_len(ng), np), , drop = FALSE],
    state$pair_info[rep(seq_len(np), each = ng), , drop = FALSE],
    replications = state$n,
    MSE_a = state$loss_mean[a], MSE_b = state$loss_mean[b],
    squared_loss_variance_a = variance[a], squared_loss_variance_b = variance[b],
    squared_loss_covariance = covariance,
    MSE_difference = state$pair_mean,
    squared_loss_difference_variance = variance_difference,
    MCSE_MSE_difference = sqrt(cky_mc_nonnegative(variance_difference) / state$n),
    RMSE_a = rmse_a, RMSE_b = rmse_b, RMSE_difference = rmse_b - rmse_a,
    MCSE_RMSE_difference = sqrt(v_rmse / state$n), row.names = NULL)
}

cky_mc_write <- function(state, output_dir, prefix, p) {
  paths <- file.path(output_dir, sprintf(c("%s_mcse-p%d.csv", "%s_paired_mcse-p%d.csv",
                                          "%s_mc_moments-p%d.rds"), prefix, p))
  write.csv(cky_mc_pointwise(state), paths[1L], row.names = FALSE)
  write.csv(cky_mc_pairs(state), paths[2L], row.names = FALSE)
  # Compact sufficient statistics permit independent checks and later merging.
  saveRDS(as.list.environment(state, all.names = TRUE), paths[3L], compress = FALSE)
  paths
}

cky_mc_metadata <- function(root) {
  list(version = 1L, helper_md5 = unname(tools::md5sum(file.path(root, "src", "monte_carlo.R"))),
       pointwise = "Welford moments of error and squared error; delta-method RMSE MCSE",
       paired = "higher depth minus lower depth; HON minus NSS; common replicate and covariance retained",
       uncertainty = "pointwise Monte Carlo uncertainty conditional on the design, not sampling confidence intervals or simultaneous bands",
       zero_rmse = "RMSE delta-method MCSE is NA if an observed RMSE is zero; all MCSEs require at least two replications",
       accumulator = "RDS stores n, grid, indices, running means, centered sums of squares and paired comoments; no replicate-level predictions")
}

cky_fold_rng <- function(seed) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    stop("Initialize the data-generation RNG before the honest-fold stream")
  }
  original <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(assign(".Random.seed", original, envir = .GlobalEnv))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  stream <- new.env(parent = emptyenv())
  stream$state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  stream$seed <- seed
  stream
}

cky_random_folds <- function(n, training_n, stream) {
  stopifnot(training_n > 0L, training_n < n)
  original <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(assign(".Random.seed", original, envir = .GlobalEnv))
  assign(".Random.seed", stream$state, envir = .GlobalEnv)
  train <- sort.int(sample.int(n, training_n, replace = FALSE))
  stream$state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  is_train <- logical(n)
  is_train[train] <- TRUE
  list(train = train, estimate = which(!is_train))
}

cky_summary_init <- function(names) {
  state <- new.env(parent = emptyenv())
  state$n <- setNames(integer(length(names)), names)
  state$missing <- state$n
  state$mean <- state$m2 <- setNames(numeric(length(names)), names)
  state$minimum <- setNames(rep(Inf, length(names)), names)
  state$maximum <- setNames(rep(-Inf, length(names)), names)
  state
}

cky_summary_update <- function(state, values) {
  stopifnot(identical(names(values), names(state$n)))
  valid <- is.finite(values)
  state$missing <- state$missing + !valid
  state$n[valid] <- state$n[valid] + 1L
  delta <- values[valid] - state$mean[valid]
  state$mean[valid] <- state$mean[valid] + delta / state$n[valid]
  state$m2[valid] <- state$m2[valid] + delta * (values[valid] - state$mean[valid])
  state$minimum[valid] <- pmin(state$minimum[valid], values[valid])
  state$maximum[valid] <- pmax(state$maximum[valid], values[valid])
  invisible(state)
}

cky_summary_table <- function(state) {
  variance <- ifelse(state$n > 1L, state$m2 / (state$n - 1L), NA_real_)
  data.frame(metric = names(state$n), n_valid = unname(state$n),
             n_undefined = unname(state$missing),
             mean = ifelse(state$n > 0L, state$mean, NA_real_),
             sd = sqrt(pmax(0, variance)), MCSE_mean = sqrt(pmax(0, variance) / state$n),
             minimum = ifelse(state$n > 0L, state$minimum, NA_real_),
             maximum = ifelse(state$n > 0L, state$maximum, NA_real_), row.names = NULL)
}

cky_fold_metrics <- function(y, d, train, estimate, xi, target) {
  one <- function(rows) {
    treated <- d[rows] == 1
    nt <- sum(treated)
    nc <- length(rows) - nt
    yt <- sum(y[rows][treated])
    yc <- sum(y[rows][!treated])
    dim <- if (nt > 0L && nc > 0L) yt / nt - yc / nc else NA_real_
    c(treated = nt, control = nc, DIM = dim,
      IPW = (yt / xi - yc / (1 - xi)) / length(rows))
  }
  construction <- one(train)
  estimation <- one(estimate)
  c(setNames(construction, paste0("training_", names(construction))),
    setNames(estimation, paste0("estimation_", names(estimation))),
    estimation_DIM_minus_target = unname(estimation["DIM"] - target),
    estimation_IPW_minus_target = unname(estimation["IPW"] - target))
}
