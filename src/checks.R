# Validate metadata, legacy RMSE tables, and Monte Carlo sufficient statistics.
# expected may be a common count or named simulation_p1/p2 and jtpa_p1/p2 counts.
cky_check_outputs <- function(root, output_dir, expected = NULL) {
  designs <- c("simulation_p1", "simulation_p2", "jtpa_p1", "jtpa_p2")
  whole <- function(x, minimum = 0) is.numeric(x) && length(x) > 0L &&
    all(is.finite(x) & x >= minimum & x == floor(x))
  if (!is.null(expected)) {
    if (!whole(expected, 2) || !(length(expected) == 1L ||
        (length(expected) == 4L && !anyDuplicated(names(expected)) &&
         setequal(names(expected), designs)))) {
      stop("expected must be a count >= 2 or four counts named ", paste(designs, collapse = ", "))
    }
  }
  require_fields <- function(x, fields, label) {
    missing <- setdiff(fields, names(x))
    if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
  require_scalars <- function(x, fields, label) {
    require_fields(x, fields, label)
    valid <- vapply(x[fields], function(value) is.atomic(value) && length(value) == 1L &&
      !is.na(value) && (!is.numeric(value) || is.finite(value)), logical(1))
    if (!all(valid)) stop(label, " has invalid scalar fields: ", paste(fields[!valid], collapse = ", "))
  }
  close <- function(a, b, label, scale = pmax(abs(a), abs(b))) {
    if (length(a) != length(b) || !identical(is.na(a), is.na(b)) ||
        any(!is.finite(a[!is.na(a)])) || any(!is.finite(b[!is.na(b)])) ||
        any(abs(a - b) > 2e-8 * pmax(1, scale), na.rm = TRUE)) {
      stop(label, " mismatch", call. = FALSE)
    }
  }
  read_table <- function(name, fields) {
    x <- read.csv(file.path(output_dir, name), check.names = FALSE)
    require_fields(x, fields, name)
    x
  }
  check_count <- function(count, design) {
    if (length(count) != 1L || !whole(count, 2)) stop(design, " needs at least two replications")
    if (!is.null(expected)) {
      wanted <- if (length(expected) == 1L) expected else expected[[design]]
      if (count != wanted) stop(design, " replication count mismatch")
    }
    count
  }
  methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE", "HON-DIM", "HON-IPW", "HON-SSE")
  depths <- 1:5
  counts <- setNames(integer(4L), designs)
  check_provenance <- function(meta, prefix, files) {
    require_fields(meta, c("status", "implementation", "script_md5", "methods",
      "monte_carlo_source_md5", "monte_carlo", "moment_file", "output_files"), prefix)
    require_scalars(meta, c("status", "script_md5", "moment_file"), prefix)
    require_scalars(meta$monte_carlo, c("version", "helper_md5", "pointwise", "paired",
      "uncertainty", "zero_rmse", "accumulator"), prefix)
    info <- meta$implementation
    require_fields(info, c("backend", "version", "source_md5"), prefix)
    require_scalars(info, c("backend", "version"), prefix)
    stopifnot(identical(meta$status, "complete"), identical(meta$methods, methods),
              identical(info$backend, "CKY canonical causal trees"), info$version == 1L)
    source_files <- c("src/estimators.R", "src/cky_theory.cpp")
    actual <- setNames(unname(tools::md5sum(file.path(root, source_files))), source_files)
    if (!identical(info$source_md5, actual)) stop(prefix, " estimator checksum mismatch")
    helper <- "src/monte_carlo.R"
    helper_md5 <- unname(tools::md5sum(file.path(root, helper)))
    if (!identical(meta$monte_carlo_source_md5, setNames(helper_md5, helper)) ||
        !identical(meta$monte_carlo$helper_md5, helper_md5)) stop(prefix, " Monte Carlo helper checksum mismatch")
    stopifnot(meta$monte_carlo$version == 1L)
    script <- if (prefix == "simulation") "CKY_2027_PNAS--simuls.R" else "CKY_2027_PNAS--empapp.R"
    if (!identical(meta$script_md5, unname(tools::md5sum(file.path(root, script))))) {
      stop(prefix, " runner checksum mismatch")
    }
    if (!setequal(meta$output_files, files) || anyDuplicated(meta$output_files) ||
        !all(file.exists(file.path(output_dir, files)))) stop(prefix, " output-file manifest mismatch")
  }
  check_diagnostics <- function(name, count) {
    fields <- c("leaves", "empty_est_leaves", "missing_arm_est_leaves",
      "zero_convention_leaves", "singleton_train_leaves", "minimum_train_size")
    z <- read_table(name, c("method", "K", "replications", fields))
    stopifnot(nrow(z) == 30L, !anyDuplicated(z[c("method", "K")]),
      setequal(z$method, methods), setequal(z$K, depths), all(z$replications == count),
      whole(as.numeric(as.matrix(z[fields]))), all(z$leaves >= count),
      all(z$leaves <= count * 2^z$K), all(z$minimum_train_size >= count),
      all(z$empty_est_leaves <= z$missing_arm_est_leaves),
      all(z$missing_arm_est_leaves <= z$leaves),
      all(z$zero_convention_leaves <= z$leaves), all(z$singleton_train_leaves <= z$leaves))
    ipw <- grepl("IPW$", z$method)
    stopifnot(all(z$zero_convention_leaves[ipw] == z$empty_est_leaves[ipw]),
      all(z$zero_convention_leaves[!ipw] == z$missing_arm_est_leaves[!ipw]),
      all(z$singleton_train_leaves[!ipw] == 0),
      all(z$zero_convention_leaves[grepl("^NSS", z$method)] == 0))
  }
  check_mc <- function(prefix, p, meta, count, grid, legacy) {
    ng <- nrow(grid)
    coords <- names(grid)
    state <- readRDS(file.path(output_dir, meta$moment_file))
    require_fields(state, c("n", "grid", "depths", "methods", "shape", "error_mean", "error_m2",
      "loss_mean", "loss_m2", "pair_info", "pair_a", "pair_b", "pair_mean", "pair_m2", "pair_comoment"), "moments")
    require_scalars(state, "n", "moments")
    stopifnot(state$n == count, identical(state$depths, depths), identical(state$methods, methods),
      identical(as.integer(state$shape), c(ng, 5L, 6L)), identical(names(state$grid), coords))
    close(as.numeric(as.matrix(state$grid)), as.numeric(as.matrix(grid)), "moment grid")
    cells <- ng * 30L
    for (field in c("error_mean", "error_m2", "loss_mean", "loss_m2")) {
      stopifnot(is.numeric(state[[field]]), length(state[[field]]) == cells, all(is.finite(state[[field]])))
    }
    stopifnot(all(state$error_m2 >= 0), all(state$loss_m2 >= 0), all(state$loss_mean >= 0))
    fields <- c("mean_error", "MCSE_mean_error", "error_variance", "MSE", "squared_loss_variance",
      "RMSE", "MCSE_MSE", "MCSE_RMSE", "relative_MCSE_RMSE")
    point <- read_table(sprintf("%s_mcse-p%d.csv", prefix, p), c(coords, "method", "K", "replications", fields))
    stopifnot(nrow(point) == cells, all(point$replications == count),
      identical(point$method, rep(methods, each = ng * 5L)),
      identical(point$K, rep(rep(depths, each = ng), 6L)),
      !anyDuplicated(point[c(coords, "method", "K")]))
    close(as.numeric(as.matrix(point[coords])), as.numeric(as.matrix(grid[rep(seq_len(ng), 30L), , drop = FALSE])), "pointwise grid")
    for (field in setdiff(fields, c("MCSE_RMSE", "relative_MCSE_RMSE"))) {
      stopifnot(is.numeric(point[[field]]), all(is.finite(point[[field]])))
    }
    stopifnot(all(as.matrix(point[setdiff(fields, "mean_error")]) >= 0, na.rm = TRUE))
    close(point$mean_error, state$error_mean, "mean error")
    close(point$error_variance, state$error_m2 / (count - 1), "error variance")
    close(point$MCSE_mean_error, sqrt(point$error_variance / count), "mean-error MCSE")
    close(point$MSE, state$loss_mean, "MSE")
    close(point$MSE, (count - 1) / count * point$error_variance + point$mean_error^2, "error second moment")
    close(point$squared_loss_variance, state$loss_m2 / (count - 1), "squared-loss variance")
    close(point$RMSE, sqrt(point$MSE), "RMSE")
    close(point$RMSE, legacy, "legacy/pointwise RMSE")
    close(point$MCSE_MSE, sqrt(point$squared_loss_variance / count), "MSE MCSE")
    delta <- ifelse(point$RMSE > 0, point$MCSE_MSE / (2 * point$RMSE), NA_real_)
    close(point$MCSE_RMSE, delta, "RMSE MCSE")
    close(point$relative_MCSE_RMSE, delta / point$RMSE, "relative RMSE MCSE")

    # Define the required comparisons independently of the saved pair indices.
    pairs <- do.call(rbind, lapply(methods, function(method) data.frame(
      comparison = "adjacent_depth", method_a = method, depth_a = 1:4,
      method_b = method, depth_b = 2:5)))
    pairs <- rbind(pairs, do.call(rbind, lapply(c("IPW", "DIM", "SSE"), function(rule) data.frame(
      comparison = "HON_minus_NSS", method_a = paste0("NSS-", rule), depth_a = depths,
      method_b = paste0("HON-", rule), depth_b = depths))))
    stopifnot(isTRUE(all.equal(state$pair_info, pairs, check.attributes = FALSE)))
    pair_rows <- pairs[rep(seq_len(nrow(pairs)), each = ng), , drop = FALSE]
    index <- function(method, depth) rep(seq_len(ng), nrow(pairs)) +
      ng * (depth - 1L) + ng * 5L * (match(method, methods) - 1L)
    a <- index(pair_rows$method_a, pair_rows$depth_a)
    b <- index(pair_rows$method_b, pair_rows$depth_b)
    stopifnot(identical(as.integer(state$pair_a), as.integer(a)), identical(as.integer(state$pair_b), as.integer(b)))
    paired_fields <- c("MSE_a", "MSE_b", "squared_loss_variance_a", "squared_loss_variance_b",
      "squared_loss_covariance", "MSE_difference", "squared_loss_difference_variance",
      "MCSE_MSE_difference", "RMSE_a", "RMSE_b", "RMSE_difference", "MCSE_RMSE_difference")
    pair <- read_table(sprintf("%s_paired_mcse-p%d.csv", prefix, p), c(coords, names(pairs), "replications", paired_fields))
    stopifnot(nrow(pair) == ng * nrow(pairs), all(pair$replications == count),
      isTRUE(all.equal(pair[names(pairs)], pair_rows, check.attributes = FALSE)))
    close(as.numeric(as.matrix(pair[coords])), as.numeric(as.matrix(grid[rep(seq_len(ng), nrow(pairs)), , drop = FALSE])), "paired grid")
    for (field in c("pair_mean", "pair_m2", "pair_comoment")) {
      stopifnot(length(state[[field]]) == nrow(pair), all(is.finite(state[[field]])))
    }
    stopifnot(all(state$pair_m2 >= 0))
    for (field in setdiff(paired_fields, "MCSE_RMSE_difference")) stopifnot(all(is.finite(pair[[field]])))
    va <- point$squared_loss_variance[a]
    vb <- point$squared_loss_variance[b]
    cv <- pair$squared_loss_covariance
    close(pair$MSE_a, point$MSE[a], "paired MSE a")
    close(pair$MSE_b, point$MSE[b], "paired MSE b")
    close(pair$squared_loss_variance_a, va, "paired variance a")
    close(pair$squared_loss_variance_b, vb, "paired variance b")
    close(cv, state$pair_comoment / (count - 1), "paired covariance")
    bound <- sqrt(va) * sqrt(vb)
    stopifnot(all(abs(cv) <= bound + 2e-8 * pmax(1, bound)))
    close(pair$MSE_difference, state$pair_mean, "paired MSE mean")
    close(pair$MSE_difference, point$MSE[b] - point$MSE[a], "paired MSE difference", pmax(point$MSE[a], point$MSE[b]))
    close(pair$squared_loss_difference_variance, state$pair_m2 / (count - 1), "paired difference variance")
    close(pair$squared_loss_difference_variance, va + vb - 2 * cv, "paired covariance identity", va + vb + 2 * abs(cv))
    stopifnot(all(pair$squared_loss_difference_variance >= 0))
    close(pair$MCSE_MSE_difference, sqrt(pair$squared_loss_difference_variance / count), "paired MSE MCSE")
    close(pair$RMSE_a, point$RMSE[a], "paired RMSE a")
    close(pair$RMSE_b, point$RMSE[b], "paired RMSE b")
    close(pair$RMSE_difference, point$RMSE[b] - point$RMSE[a], "paired RMSE difference")
    # Use unrounded moments for the delta expansion: cancellation can amplify
    # the CSV's decimal rounding when paired predictions nearly coincide.
    ma <- state$loss_mean[a]
    mb <- state$loss_mean[b]
    positive <- ma > 0 & mb > 0
    delta_var <- rep(NA_real_, nrow(pair))
    scale <- state$loss_m2[a[positive]] / (count - 1) / (4 * ma[positive]) +
      state$loss_m2[b[positive]] / (count - 1) / (4 * mb[positive])
    value <- scale - state$pair_comoment[positive] / (count - 1) / (2 * sqrt(ma[positive]) * sqrt(mb[positive]))
    stopifnot(all(value >= -2e-8 * pmax(1, scale)))
    delta_var[positive] <- pmax(0, value) / count
    close(pair$MCSE_RMSE_difference, sqrt(delta_var), "paired RMSE MCSE")
  }
  check_folds <- function(meta, count, raw) {
    stopifnot(identical(meta$fold_rule, paste0("fresh uniform subset of floor(n/2) rows for construction in every replication; ",
        "complementary rows for estimation; no outcome- or treatment-based rejection")),
      identical(meta$fold_rng_kind, c("L'Ecuyer-CMRG", "Inversion", "Rejection")),
      isTRUE(meta$fold_rng_separate_from_covariate_permutations), meta$fold_seed == 12000L + meta$p,
      meta$training_n == floor(nrow(raw) / 2), meta$estimation_n == nrow(raw) - meta$training_n)
    metrics <- c(paste0("training_", c("treated", "control", "DIM", "IPW")),
      paste0("estimation_", c("treated", "control", "DIM", "IPW")),
      "estimation_DIM_minus_target", "estimation_IPW_minus_target")
    z <- read_table(meta$fold_diagnostics_file, c("metric", "n_valid", "n_undefined", "mean", "sd", "MCSE_mean", "minimum", "maximum"))
    stopifnot(nrow(z) == length(metrics), !anyDuplicated(z$metric), setequal(z$metric, metrics),
      whole(z$n_valid), whole(z$n_undefined), all(z$n_valid + z$n_undefined == count))
    for (j in seq_len(nrow(z))) {
      values <- unlist(z[j, c("mean", "minimum", "maximum")], use.names = FALSE)
      if (z$n_valid[j] == 0L) stopifnot(all(is.na(values))) else {
        stopifnot(all(is.finite(values)), z$minimum[j] <= z$mean[j], z$mean[j] <= z$maximum[j])
      }
      if (z$n_valid[j] < 2L) stopifnot(is.na(z$sd[j]), is.na(z$MCSE_mean[j])) else {
        stopifnot(is.finite(z$sd[j]), z$sd[j] >= 0)
        close(z$MCSE_mean[j], z$sd[j] / sqrt(z$n_valid[j]), "fold-summary MCSE")
      }
    }
    get <- function(metric, field = "mean") z[[field]][match(metric, z$metric)]
    for (fold in c("training", "estimation")) {
      size <- meta[[paste0(fold, "_n")]]
      for (arm in c("treated", "control")) {
        metric <- paste0(fold, "_", arm)
        stopifnot(get(metric, "n_valid") == count, get(metric, "minimum") >= 0,
          get(metric, "maximum") <= size, whole(c(get(metric, "minimum"), get(metric, "maximum"))))
      }
      close(get(paste0(fold, "_treated")) + get(paste0(fold, "_control")), size, "fold count sum")
      close(get(paste0(fold, "_treated"), "sd"), get(paste0(fold, "_control"), "sd"), "fold count SD")
      stopifnot(get(paste0(fold, "_IPW"), "n_valid") == count)
    }
    for (arm in c("treated", "control")) {
      total <- meta[[paste0("n_", arm)]]
      train <- paste0("training_", arm)
      est <- paste0("estimation_", arm)
      close(get(train) + get(est), total, "complementary arm count")
      close(get(train, "sd"), get(est, "sd"), "complementary count SD")
      close(get(train, "minimum") + get(est, "maximum"), total, "complementary count extrema")
      close(get(train, "maximum") + get(est, "minimum"), total, "complementary count extrema")
    }
    close(meta$training_n * get("training_IPW") + meta$estimation_n * get("estimation_IPW"),
          nrow(raw) * meta$pooled_IPW, "complementary IPW means")
    close(meta$training_n * get("training_IPW", "sd"), meta$estimation_n * get("estimation_IPW", "sd"), "complementary IPW SD")
    for (estimator in c("DIM", "IPW")) {
      basic <- paste0("estimation_", estimator)
      difference <- paste0(basic, "_minus_target")
      stopifnot(get(basic, "n_valid") == get(difference, "n_valid"))
      for (field in c("mean", "minimum", "maximum")) close(get(difference, field), get(basic, field) - meta$tau, "fold target offset")
      close(get(difference, "sd"), get(basic, "sd"), "fold target-offset SD")
    }
  }
  for (prefix in c("simulation", "jtpa")) for (p in 1:2) {
    design <- paste0(prefix, "_p", p)
    meta_name <- if (prefix == "simulation") sprintf("simulation_metadata-p%d.R", p) else sprintf("jtpa_p%d_metadata.R", p)
    meta <- dget(file.path(output_dir, meta_name))
    require_fields(meta, c("moment_file", "rng_kind"), meta_name)
    if (prefix == "simulation") {
      require_fields(meta, c("config", "design"), meta_name)
      require_scalars(meta, "counts_per_grid_point_and_depth", meta_name)
      require_scalars(meta$config, c("p", "n", "seed", "m_reps"), meta_name)
      require_fields(meta$config, c("K_list", "grid"), meta_name)
      require_scalars(meta$design, c("tau", "treatment_probability"), meta_name)
    } else {
      require_scalars(meta, c("p", "seed", "m_reps", "n_valid", "n_input", "n_dropped_missing", "n_complete",
        "n_treated", "n_control", "training_n", "estimation_n", "input_file", "input_md5", "tau",
        "pooled_DIM", "pooled_IPW", "treatment_probability", "grid_points", "fold_rule", "fold_seed",
        "fold_rng_separate_from_covariate_permutations", "fold_diagnostics_file"), meta_name)
      require_fields(meta, c("depths", "fold_rng_kind"), meta_name)
    }
    count <- check_count(if (prefix == "simulation") meta$config$m_reps else meta$m_reps, design)
    counts[[design]] <- count
    coords <- paste0("X", seq_len(p))
    mc_files <- sprintf(c("%s_mcse-p%d.csv", "%s_paired_mcse-p%d.csv", "%s_mc_moments-p%d.rds"), prefix, p)
    stopifnot(identical(meta$moment_file, mc_files[3L]),
      identical(meta$rng_kind, c("Mersenne-Twister", "Inversion", "Rejection")))
    diagnostic <- sprintf("%s_diagnostics-p%d.csv", prefix, p)
    if (prefix == "simulation") {
      stopifnot(meta$config$p == p, meta$config$n == 1000L, meta$config$seed == 123L,
        identical(meta$config$K_list, depths), meta$counts_per_grid_point_and_depth == count,
        meta$design$tau == 1, meta$design$treatment_probability == .5)
      grid <- if (p == 1L) data.frame(X1 = seq(0, 1, length.out = 201)) else
        expand.grid(X1 = seq(0, 1, length.out = 51), X2 = seq(0, 1, length.out = 51))
      close(as.numeric(as.matrix(meta$config$grid)), as.numeric(as.matrix(grid)), "simulation metadata grid")
      legacy <- numeric()
      old_files <- sprintf("rmse_%s-p%d.csv", gsub("-", "_", methods), p)
      for (j in seq_along(methods)) {
        x <- read_table(old_files[j], c("p", "method", "n", "m_reps", "tau", "seed", coords, paste0("K", depths), paste0("N_K", depths)))
        stopifnot(nrow(x) == nrow(grid), all(x$p == p), all(x$method == methods[j]),
          all(x$n == 1000L), all(x$m_reps == count), all(x$tau == 1), all(x$seed == 123L),
          all(as.matrix(x[paste0("N_K", depths)]) == count))
        close(as.numeric(as.matrix(x[coords])), as.numeric(as.matrix(grid)), "simulation CSV grid")
        legacy <- c(legacy, as.numeric(as.matrix(x[paste0("K", depths)])))
      }
      check_provenance(meta, prefix, c(old_files, diagnostic, mc_files))
    } else {
      input_name <- if (p == 1L) "jtpa_prevearn_0_6000_subsample.csv" else "jtpa_prevearn_0_5000_subsample_p2.csv"
      input <- file.path(root, "jtpa", input_name)
      raw <- read.csv(input)
      n_input <- nrow(raw)
      raw <- raw[complete.cases(raw[c("y", "d", coords)]), , drop = FALSE]
      stopifnot(meta$p == p, meta$seed == 12L, meta$n_valid == count, identical(meta$depths, depths),
        meta$n_input == n_input, meta$n_dropped_missing == n_input - nrow(raw), meta$n_complete == nrow(raw),
        meta$n_treated == sum(raw$d == 1), meta$n_control == sum(raw$d == 0),
        identical(meta$input_file, input_name), identical(meta$input_md5, unname(tools::md5sum(input))))
      close(meta$treatment_probability, 2/3, "treatment probability")
      tau <- mean(raw$y[raw$d == 1]) - mean(raw$y[raw$d == 0])
      close(meta$tau, tau, "empirical target")
      close(meta$pooled_DIM, tau, "pooled DIM")
      close(meta$pooled_IPW, mean(raw$y * (raw$d - 2/3) / ((2/3) * (1/3))), "pooled IPW")
      if (p == 1L) {
        values <- sort(unique(raw$X1))
        if (length(values) > 201L) values <- values[unique(round(seq(1, length(values), length.out = 201)))]
        grid <- data.frame(X1 = values)
      } else grid <- expand.grid(
        X1 = seq(quantile(raw$X1, .01), quantile(raw$X1, .99), length.out = 41),
        X2 = seq(quantile(raw$X2, .01), quantile(raw$X2, .99), length.out = 31))
      stopifnot(meta$grid_points == nrow(grid))
      old_file <- if (p == 1L) "rmse_K1to5.csv" else "rmse_K1to5_p2.csv"
      x <- read_table(old_file, c(coords, "method", "K", "RMSE", "n", "m_reps", "tau", "seed", "n_valid"))
      stopifnot(nrow(x) == 30L * nrow(grid), all(x$n == nrow(raw)), all(x$m_reps == count),
        all(x$n_valid == count), all(x$seed == 12L),
        identical(x$method, rep(methods, each = 5L * nrow(grid))),
        identical(x$K, rep(rep(depths, each = nrow(grid)), 6L)))
      close(x$tau, rep(tau, nrow(x)), "CSV empirical target")
      close(as.numeric(as.matrix(x[coords])), as.numeric(as.matrix(grid[rep(seq_len(nrow(grid)), 30L), , drop = FALSE])), "empirical CSV grid")
      legacy <- x$RMSE
      stopifnot(identical(meta$fold_diagnostics_file, sprintf("jtpa_fold_diagnostics-p%d.csv", p)))
      check_provenance(meta, prefix, c(old_file, diagnostic, mc_files, meta$fold_diagnostics_file))
      check_folds(meta, count, raw)
    }
    check_diagnostics(diagnostic, count)
    check_mc(prefix, p, meta, count, grid, legacy)
  }
  cat("Validated all four experiments, six methods, five depths, and Monte Carlo statistics.\n",
      paste(names(counts), counts, sep = "=", collapse = ", "), " replications.\n", sep = "")
  invisible(TRUE)
}

# Independent oracle: enumerate candidate partitions; SSE uses OLS residual sums.
cky_test_engine <- function(root) {
  cky_load_engine(root)
  checks <- 0L
  assert <- function(ok, description) {
    if (!isTRUE(ok)) stop(description, call. = FALSE)
    checks <<- checks + 1L
  }
  close <- function(a, b, description) {
    assert(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)), description)
  }
  oracle <- function(x, y, d, rule, xi, K, ex, ey, ed, grid) {
    leaf <- function(rows) {
      if (rule == "IPW") return(if (length(rows)) mean(ey[rows] * (ed[rows] - xi) / (xi * (1-xi))) else 0)
      t <- rows[ed[rows] == 1]; c <- rows[ed[rows] == 0]
      if (!length(t) || !length(c)) 0 else mean(ey[t]) - mean(ey[c])
    }
    fitmean <- function(rows) {
      if (rule == "IPW") mean(y[rows] * (d[rows] - xi) / (xi*(1-xi))) else
        mean(y[rows[d[rows] == 1]]) - mean(y[rows[d[rows] == 0]])
    }
    rss <- function(rows) {
      # Direct residual fit of y on intercept and treatment, independently in a child.
      residual <- lm.fit(cbind(1, d[rows]), y[rows])$residuals
      sum(residual^2)
    }
    grow <- function(rows, est, points, depth) {
      best <- -Inf; chosen <- NULL
      if (depth < K && length(rows) > 1L) {
        for (j in seq_len(ncol(x))) {
          cuts <- head(sort(unique(x[rows,j])), -1L)
          for (cut in cuts) {
            L <- rows[x[rows,j] <= cut]; R <- rows[x[rows,j] > cut]
            if (rule != "IPW" && (!all(c(0,1) %in% d[L]) || !all(c(0,1) %in% d[R]))) next
            score <- if (rule == "SSE") rss(rows)-rss(L)-rss(R) else
              length(L)*length(R)/length(rows)*(fitmean(L)-fitmean(R))^2
            if (!is.finite(best) || score > best + 64*.Machine$double.eps*max(abs(best),abs(score))) {best <- score; chosen <- list(j=j,cut=cut,L=L,R=R)}
          }
        }
      }
      if (is.null(chosen)) return(rep(leaf(est), length(points)))
      left <- grid[points,chosen$j] <= chosen$cut
      est_left <- ex[est,chosen$j] <= chosen$cut
      answer <- numeric(length(points))
      answer[left] <- grow(chosen$L, est[est_left], points[left], depth+1L)
      answer[!left] <- grow(chosen$R, est[!est_left], points[!left], depth+1L)
      answer
    }
    grow(seq_len(nrow(x)), seq_len(nrow(ex)), seq_len(nrow(grid)), 0L)
  }

  # IPW is a transformed-outcome mean even when local arm proportions differ.
  x <- matrix(1:4, ncol=1); y <- c(1,1,1,0); d <- c(1,1,1,0)
  ipw <- cky_fit(x,y,d,"IPW",.5,0)
  dim <- cky_fit(x,y,d,"DIM",.5,0)
  close(predict(ipw,x), rep(1.5,4), "IPW leaf mean must be 1.5")
  close(predict(dim,x), rep(1,4), "DIM leaf mean must be 1")

  # One-arm IPW construction leaves are admissible; DIM/SSE omit those splits.
  tiny_x <- matrix(1:2,ncol=1)
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(tiny_x,c(1,0),c(1,0),rule,.5,1)
    assert(nrow(tr$nodes) == if(rule == "IPW") 3L else 1L, paste(rule,"admissibility"))
  }

  # Empty and missing-arm honest estimation leaves return the specified values.
  tr <- cky_fit(tiny_x,c(1,0),c(1,0),"IPW",.5,1,
                est_x=matrix(1,ncol=1),est_y=3,est_d=1)
  close(predict(tr,tiny_x),c(6,0),"IPW uses one-arm mean and zero only when empty")
  for (rule in c("DIM","SSE")) {
    tx <- matrix(1:4,ncol=1); ty <- c(2,0,8,1); td <- c(1,0,1,0)
    tr <- cky_fit(tx,ty,td,rule,.5,1,est_x=matrix(c(1,3,4),ncol=1),
                  est_y=c(20,10,2),est_d=c(1,1,0))
    close(predict(tr,tx),c(0,0,8,8),paste(rule,"must not inherit the parent estimate"))
  }
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(x,y,d,rule,.5,3,est_x=matrix(numeric(),ncol=1),est_y=numeric(),est_d=numeric())
    close(predict(tr,x),rep(0,4),paste(rule,"entire estimation sample empty"))
  }

  # Observed-value thresholds, <= convention, coordinate ties, zero-gain splits.
  tr <- cky_fit(cbind(c(0,2),c(0,2)),c(3,0),c(1,0),"IPW",.5,1)
  assert(tr$nodes[1,"feature"]==1 && tr$nodes[1,"cut"]==0,"First coordinate and observed-value cutoff")
  close(predict(tr,cbind(c(0,1,2),c(0,1,2))),c(6,0,0),"Prediction between observations goes right of observed cutoff")
  for (rule in c("IPW","DIM","SSE")) {
    tr <- cky_fit(matrix(1:8,ncol=1),rep(0,8),rep(0:1,4),rule,.5,2)
    assert(nrow(tr$nodes)>1 && tr$nodes[1,"gain"]==0,paste(rule,"no positive-gain stopping restriction"))
    same_x <- cky_fit(matrix(1,8,2),1:8,rep(0:1,4),rule,.5,3)
    assert(nrow(same_x$nodes)==1,paste(rule,"no splitting equal covariates"))
  }

  # Hundreds of small comparisons against independently enumerated direct formulas.
  set.seed(97431)
  for (trial in 1:70) {
    n <- sample(8:30,1); p <- sample(1:3,1); ne <- sample(0:25,1)
    x <- matrix(if (trial %% 2) runif(n*p) else sample(0:9,n*p,TRUE),n,p)
    y <- rnorm(n,mean=2); d <- rbinom(n,1,.4)
    ex <- matrix(if (trial %% 2) runif(ne*p) else sample(0:9,ne*p,TRUE),ne,p)
    ey <- rnorm(ne); ed <- rbinom(ne,1,.6)
    grid <- rbind(x,ex,matrix(rep(c(-1,11),p),2,p))
    for (rule in c("IPW","DIM","SSE")) {
      tr <- cky_fit(x,y,d,rule,.4,3,est_x=ex,est_y=ey,est_d=ed)
      pred <- predict(tr,grid,depths=0:3)
      for (K in 0:3) {
        ref <- oracle(x,y,d,rule,.4,K,ex,ey,ed,grid)
        close(pred[,K+1],ref,paste("Oracle mismatch",trial,rule,K))
        separate <- cky_fit(x,y,d,rule,.4,K,est_x=ex,est_y=ey,est_d=ed)
        close(pred[,K+1],predict(separate,grid),paste("Depth truncation mismatch",trial,rule,K))
      }
      changed <- cky_fit(x,y,d,rule,.4,3,est_x=ex,est_y=ey+100,est_d=1-ed)
      columns <- c("left","right","feature","cut","depth","n_train","n1_train","gain","n0_train")
      assert(identical(tr$nodes[,columns],changed$nodes[,columns]),paste(rule,"honest outcomes leaked into construction"))
    }
  }

  # Fixed DIM objective fixture: the maximizing cutoff is 6.
  fixture_y <- c(-2,4,0,-3,-2,0,3,0,0,4,1,0)
  fixture_d <- c(0,0,0,0,1,1,0,1,1,0,1,1)
  fixture_x <- matrix(1:12,ncol=1)
  fixture <- cky_fit(fixture_x,fixture_y,fixture_d,"DIM",.5,1)
  assert(fixture$nodes[1,"cut"]==6,"DIM must select the maximizing child-contrast cutoff")
  for (scale in c(1e-9,1e9)) {
    scaled <- cky_fit(fixture_x,fixture_y*scale,fixture_d,"DIM",.5,1)
    assert(scaled$nodes[1,"cut"]==6,"Numerical tie rule must preserve outcome scaling")
  }

  # No stochastic fits: calls must not consume the simulation RNG stream.
  before <- .Random.seed
  invisible(cky_fit(x,y,d,"IPW",.4,3))
  assert(identical(before,.Random.seed),"Tree fitting consumed RNG")
  assert(inherits(try(cky_fit(x,y,d,"IPW",0,1),silent=TRUE),"try-error"),"Reject invalid propensity")
  assert(inherits(try(cky_fit(x,y,d,"DIM",.4,-1),silent=TRUE),"try-error"),"Reject invalid depth")
  cat("Passed", checks,"checks against direct formulas, edge cases, honest isolation, and depth equivalence.\n")
  invisible(TRUE)
}

# Batch references intentionally use base mean/var/cov on saved fixture draws.
cky_test_monte_carlo <- function(root) {
  source(file.path(root, "src", "monte_carlo.R"), local = environment())
  checks <- 0L
  assert <- function(ok, label) {
    if (!isTRUE(ok)) stop(label, call. = FALSE)
    checks <<- checks + 1L
  }
  close <- function(a, b, label) {
    assert(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)), label)
  }
  methods <- c("NSS-DIM", "NSS-IPW", "NSS-SSE", "HON-DIM", "HON-IPW", "HON-SSE")
  grid <- data.frame(X1 = c(.1, .9))
  depths <- 1:3
  m <- 7L
  # Vary signs, means, scales, and correlations across depth/method; one entire
  # grid point is zero so zero-boundary pairs occur alongside ordinary pairs.
  errors <- array(0, c(2L, 3L, 6L, m), dimnames = list(NULL, as.character(depths), methods, NULL))
  for (j in seq_along(methods)) for (k in depths) {
    errors[2L, k, j, ] <- (seq_len(m) - 3L) * (j + k) / 5 + (-1)^seq_len(m) * (j - k) / 7
  }
  # Exact proportional losses provide a strongly correlated paired fixture.
  errors[2L, 1L, "HON-DIM", ] <- 2 * errors[2L, 1L, "NSS-DIM", ]
  check_fixture <- function(draws, label) {
    state <- cky_mc_init(grid, depths, methods)
    target <- 1.25
    for (r in seq_len(dim(draws)[4L])) cky_mc_update(state, draws[, , , r] + target, target)
    point <- cky_mc_pointwise(state)
    loss <- draws^2
    batch_mean <- as.numeric(apply(draws, 1:3, mean))
    batch_var <- as.numeric(apply(draws, 1:3, var))
    loss_mean <- as.numeric(apply(loss, 1:3, mean))
    loss_var <- as.numeric(apply(loss, 1:3, var))
    nr <- dim(draws)[4L]
    close(point$mean_error, batch_mean, paste(label, "mean error"))
    close(point$error_variance, batch_var, paste(label, "error variance"))
    close(point$MCSE_mean_error, sqrt(batch_var / nr), paste(label, "mean-error MCSE"))
    close(point$MSE, loss_mean, paste(label, "MSE"))
    close(point$squared_loss_variance, loss_var, paste(label, "squared-loss variance"))
    close(point$RMSE, sqrt(loss_mean), paste(label, "RMSE"))
    close(point$MCSE_MSE, sqrt(loss_var / nr), paste(label, "MSE MCSE"))
    delta <- ifelse(loss_mean > 0, sqrt(loss_var / nr) / (2 * sqrt(loss_mean)), NA_real_)
    close(point$MCSE_RMSE, delta, paste(label, "RMSE MCSE"))
    close(point$relative_MCSE_RMSE, delta / sqrt(loss_mean), paste(label, "relative MCSE"))
    pair <- cky_mc_pairs(state)
    assert(nrow(pair) == 2L * (6L * 2L + 3L * 3L), paste(label, "pair coverage"))
    for (i in seq_len(nrow(pair))) {
      g <- match(pair$X1[i], grid$X1)
      a <- loss[g, as.character(pair$depth_a[i]), pair$method_a[i], ]
      b <- loss[g, as.character(pair$depth_b[i]), pair$method_b[i], ]
      diff <- b - a
      close(pair$MSE_a[i], mean(a), paste(label, "paired mean a"))
      close(pair$MSE_b[i], mean(b), paste(label, "paired mean b"))
      close(pair$squared_loss_variance_a[i], var(a), paste(label, "paired variance a"))
      close(pair$squared_loss_variance_b[i], var(b), paste(label, "paired variance b"))
      close(pair$squared_loss_covariance[i], cov(a, b), paste(label, "paired covariance"))
      close(pair$MSE_difference[i], mean(diff), paste(label, "paired mean difference"))
      close(pair$squared_loss_difference_variance[i], var(diff), paste(label, "paired difference variance"))
      close(pair$MCSE_MSE_difference[i], sd(diff) / sqrt(nr), paste(label, "paired MSE MCSE"))
      close(pair$RMSE_difference[i], sqrt(mean(b)) - sqrt(mean(a)), paste(label, "paired RMSE difference"))
      # Direct batch influence draws avoid reusing the covariance expansion.
      direct <- if (mean(a) > 0 && mean(b) > 0) {
        influence <- (b - mean(b)) / (2 * sqrt(mean(b))) - (a - mean(a)) / (2 * sqrt(mean(a)))
        sd(influence) / sqrt(nr)
      } else NA_real_
      close(pair$MCSE_RMSE_difference[i], direct, paste(label, "paired RMSE MCSE"))
    }
    invisible(state)
  }
  mixed <- check_fixture(errors, "mixed/correlated")
  check_fixture(errors * 0, "all zero")
  single <- cky_mc_init(grid, depths, methods)
  cky_mc_update(single, errors[, , , 1L], 0)
  one <- cky_mc_pointwise(single)
  assert(all(is.na(one$MCSE_mean_error)) && all(is.na(one$MCSE_MSE)) && all(is.na(one$MCSE_RMSE)),
    "one replication leaves MCSE undefined")
  assert(all(is.na(cky_mc_pairs(single)$MCSE_MSE_difference)), "one replication leaves paired MCSE undefined")

  # Serialization must retain the accumulator rather than rounded CSV values.
  scratch <- tempfile("cky-moments-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  saved <- cky_mc_write(mixed, scratch, "fixture", 1L)
  restored <- list2env(readRDS(saved[3L]), parent = emptyenv())
  assert(identical(cky_mc_pointwise(mixed), cky_mc_pointwise(restored)), "moment RDS pointwise round trip")
  assert(identical(cky_mc_pairs(mixed), cky_mc_pairs(restored)), "moment RDS paired round trip")

  values <- cbind(a = c(1, 4, NA, -2, 8), b = c(0, 0, 0, 0, 0), c = c(NA, Inf, 3, NA, NA))
  summary <- cky_summary_init(colnames(values))
  for (j in seq_len(nrow(values))) cky_summary_update(summary, values[j, ])
  table <- cky_summary_table(summary)
  for (j in seq_len(ncol(values))) {
    direct <- values[is.finite(values[, j]), j]
    assert(table$n_valid[j] == length(direct) && table$n_undefined[j] == nrow(values) - length(direct), "summary finite counts")
    close(table$mean[j], mean(direct), "summary mean")
    close(table$sd[j], sd(direct), "summary SD")
    close(table$MCSE_mean[j], sd(direct) / sqrt(length(direct)), "summary MCSE")
    close(table$minimum[j], min(direct), "summary minimum")
    close(table$maximum[j], max(direct), "summary maximum")
  }

  # Compare to explicit sampling on a private stream, and verify that changing
  # the number of intervening folds cannot alter the main permutation stream.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(12001L)
  direct_folds <- replicate(12L, sort(sample.int(17L, 8L, replace = FALSE)), simplify = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(12L)
  main_state <- .Random.seed
  direct_main <- replicate(12L, sample.int(17L), simplify = FALSE)
  assign(".Random.seed", main_state, envir = .GlobalEnv)
  stream <- cky_fold_rng(12001L)
  assert(identical(main_state, .Random.seed), "fold initialization preserves main RNG state")
  assert(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")), "fold initialization preserves main RNG kind")
  for (j in seq_len(12L)) {
    before <- .Random.seed
    fold <- cky_random_folds(17L, 8L, stream)
    assert(identical(before, .Random.seed), "fold draw preserves main RNG state")
    assert(identical(fold$train, direct_folds[[j]]), "fold draw matches direct private-stream sampling")
    assert(length(fold$train) == 8L && length(fold$estimate) == 9L &&
      !length(intersect(fold$train, fold$estimate)) && identical(sort(c(fold$train, fold$estimate)), 1:17), "folds partition all rows")
    assert(identical(sample.int(17L), direct_main[[j]]), "main permutation stream unchanged")
  }
  assert(length(unique(vapply(direct_folds, paste, character(1), collapse = ","))) > 1L, "folds change across replications")
  y <- c(2, 7, 4, 11)
  d <- c(0, 1, 0, 1)
  metric <- cky_fold_metrics(y, d, 1:2, 3:4, .5, 6)
  close(metric[c("training_treated", "training_control", "estimation_treated", "estimation_control")], rep(1, 4), "fold arm counts")
  close(metric[c("training_DIM", "training_IPW", "estimation_DIM", "estimation_IPW")], c(5, 5, 7, 7), "fold root contrasts")
  close(metric[c("estimation_DIM_minus_target", "estimation_IPW_minus_target")], c(1, 1), "fold target offsets")
  one_arm <- cky_fold_metrics(y, d, c(1L, 3L), c(2L, 4L), .5, 6)
  assert(is.na(one_arm[["training_DIM"]]) && is.na(one_arm[["estimation_DIM"]]) &&
    is.finite(one_arm[["training_IPW"]]) && is.finite(one_arm[["estimation_IPW"]]), "one-arm DIM diagnostics undefined, IPW defined")
  cat("Passed", checks, "independent Monte Carlo and fold checks.\n")
  invisible(checks)
}

# The manifest uses repository-relative paths, including the bundled outputs.
cky_manifest <- function(root, write = FALSE) {
  oldwd <- setwd(root)
  on.exit(setwd(oldwd), add = TRUE)
  manifest <- file.path("output", "SHA256SUMS")
  if (write) {
    files <- c("CKY_2027_PNAS--simuls.R", "CKY_2027_PNAS--empapp.R",
               "CKY_2027_PNAS--plots.R", "src/estimators.R", "src/checks.R",
               "src/cky_theory.cpp", "src/monte_carlo.R",
               "jtpa/prepare.R",
               "jtpa/jtpa_prevearn_0_6000_subsample.csv",
               "jtpa/jtpa_prevearn_0_5000_subsample_p2.csv",
               list.files("output", pattern = "\\.(csv|R|txt|rds)$", full.names = TRUE))
    files <- sort(files)
    sums <- system2("sha256sum", shQuote(files), stdout = TRUE)
    if (!is.null(attr(sums, "status")) || length(sums) != length(files)) {
      stop("Could not calculate package checksums")
    }
    writeLines(sums, manifest)
    cat("Saved", length(files), "checksums to output/SHA256SUMS\n")
  } else {
    result <- system2("sha256sum", c("--check", "--status", shQuote(manifest)))
    if (result != 0L) stop("Package checksum verification failed")
    cat("Package checksums verified.\n")
  }
  invisible(TRUE)
}

# Optional verification utility; the three main scripts remain the public workflow.
if (sys.nframe() == 0L) {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script) != 1L) stop("Run with Rscript src/checks.R [check|test]")
  root <- dirname(dirname(normalizePath(sub("^--file=", "", script))))
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[1] else "check"
  if (length(args) > 1L || !mode %in% c("check", "test")) {
    stop("Use Rscript src/checks.R [check|test]")
  }
  source(file.path(root, "src", "estimators.R"))
  if (mode == "test") {
    cky_test_engine(root)
    cky_test_monte_carlo(root)
  } else {
    output_dir <- path.expand(Sys.getenv("CKY_OUTPUT_DIR", "output"))
    if (!grepl("^/", output_dir)) output_dir <- file.path(root, output_dir)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    cky_check_outputs(root, output_dir)
    if (identical(output_dir, normalizePath(file.path(root, "output"), mustWork = TRUE))) {
      cky_manifest(root)
    } else cat("Custom output directory validated; canonical manifest not checked.\n")
  }
}
