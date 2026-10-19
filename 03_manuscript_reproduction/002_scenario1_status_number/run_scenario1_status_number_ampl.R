library(rAMPL)

# rAMPL version of run_scenario1_status_number.R.
#
# The optimization model is no longer assembled in R for the Gurobi API.
# Instead, the AMPL formulation in 01_models/001_core/OptOrch_core.mod is
# loaded through rAMPL (see 01_models/optorch_ampl.R), the data are pushed from
# R, and Gurobi is called through AMPL with the same tuned parameters.
# Output files have the same names and columns as the original script, plus
# two extra columns in summary.csv (solve_result, solve_message).

#---------------------------------------------------------------
# 0) Path ------------------------------------------------------

# Set the working directory to:
# 03_manuscript_reproduction/002_scenario1_status_number/

SCRIPT_DIR <- getwd()

DATA_DIR <- file.path(SCRIPT_DIR, "..", "001_simulated_data", "MoBPS_generated_data", "Heri_0.2")
OUT_DIR <- file.path(SCRIPT_DIR, "outputs")
MODEL_DIR <- file.path(SCRIPT_DIR, "..", "..", "01_models")

source(file.path(MODEL_DIR, "optorch_ampl.R"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Remove previous outputs if needed.
# Be careful: this removes all files and folders inside OUT_DIR, including file_description.md.
# unlink(list.files(OUT_DIR, full.names = TRUE), recursive = TRUE, force = TRUE)

#---------------------------------------------------------------
# 1) Parameter -------------------------------------------------

nc <- 500
ns_grid <- c(10, 20, 30, 40)
iter_grid <- 1:30

# Optional overrides for quick tests, e.g. in a shell:
#   OPTORCH_NS_GRID=30 OPTORCH_ITER_GRID=1 OPTORCH_TIME_LIMIT=120 Rscript run_scenario1_status_number.R
if (nzchar(Sys.getenv("OPTORCH_NS_GRID"))) {
  ns_grid <- as.integer(strsplit(Sys.getenv("OPTORCH_NS_GRID"), ",")[[1]])
}
if (nzchar(Sys.getenv("OPTORCH_ITER_GRID"))) {
  iter_grid <- as.integer(strsplit(Sys.getenv("OPTORCH_ITER_GRID"), ",")[[1]])
}
time_limit <- if (nzchar(Sys.getenv("OPTORCH_TIME_LIMIT"))) as.numeric(Sys.getenv("OPTORCH_TIME_LIMIT")) else NULL

L <- rep(0.01, nc)
U <- rep(0.15, nc)

# Tuned Gurobi parameters used for each status number (native Gurobi names,
# passed through the AMPL Gurobi driver unchanged).
tuned_params_by_ns <- list(
  "10" = list(
    MIPGap = 0.005,
    OutputFlag = 1,
    NonConvex = 2,
    Aggregate = 0,
    NormAdjust = 0,
    NumericFocus = 3,
    PreMIQCPForm = 1,
    ScaleFlag = 2
  ),
  "20" = list(
    MIPGap = 0.005,
    OutputFlag = 1,
    NonConvex = 2,
    Heuristics = 0.5,
    MIQCPMethod = 0,
    PrePasses = 1,
    Presolve = 1,
    ScaleFlag = 0
  ),
  "30" = list(
    MIPGap = 0.005,
    OutputFlag = 1,
    NonConvex = 2,
    Aggregate = 0,
    CutPasses = 10,
    Cuts = 1,
    PreMIQCPForm = 1,
    Presolve = 1,
    ScaleFlag = 0,
    ZeroHalfCuts = 0
  ),
  "40" = list(
    MIPGap = 0.005,
    OutputFlag = 1,
    NonConvex = 2,
    Heuristics = 0.001,
    NoRelHeurTime = 900,
    ScaleFlag = 2
  )
)

#---------------------------------------------------------------
# 2) AMPL session ----------------------------------------------

# AMPL is found on PATH, or in the folder given by the AMPL_HOME environment variable.
ampl <- optorch_ampl_start(solver = "gurobi")
on.exit(ampl$close(), add = TRUE)

optorch_load_model(ampl, MODEL_DIR, scenario = "status_number")

#---------------------------------------------------------------
# 3) Load Data and Run Optimization ----------------------------

for (ns in ns_grid) {

  ns_dir <- file.path(OUT_DIR, sprintf("ns_%02d", ns))
  dir.create(ns_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n===== ns = %d started =====\n", ns))
  flush.console()

  for (iter in iter_grid) {

    cat(sprintf("ns=%d | iter=%d running...\n", ns, iter))
    flush.console()

    iter_dir <- file.path(ns_dir, paste0("rep_", iter))
    dir.create(iter_dir, showWarnings = FALSE, recursive = TRUE)

    results_summary_iter <- list()
    row_id_iter <- 0

    selected_rows_iter <- list()
    sel_row_id_iter <- 0

    c_path <- list.files(
      DATA_DIR,
      pattern = paste0("^rep", iter, "_Gmatrix_tuned_all\\.csv$"),
      full.names = TRUE
    )[1]

    b_path <- list.files(
      DATA_DIR,
      pattern = paste0("^rep", iter, "_GBLUP_results\\.csv$"),
      full.names = TRUE
    )[1]

    if (is.na(c_path) || is.na(b_path)) {
      cat(sprintf("ns=%d iter=%d skipped (missing C or BV file)\n", ns, iter))
      next
    }

    # C matrix: rows and columns 102-601 in the original CSV
    dfC <- read.csv(
      c_path,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    C_block <- dfC[102:601, 102:601, drop = FALSE]
    C_full <- apply(C_block, 2, function(x) suppressWarnings(as.numeric(x)))
    C_full <- as.matrix(C_full)

    # Convert the genomic relationship matrix G to the coancestry matrix C.
    C_mat <- C_full * 0.5

    # Force PSD
    # Replace eigenvalues below 1e-6 and reconstruct the matrix
    eig <- eigen(C_mat, symmetric = TRUE)
    eig$values[eig$values < 1e-6] <- 1e-6
    C_mat <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

    # Breeding values
    dfB <- read.csv(
      b_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    raw_id <- as.character(dfB[[1]][101:600])
    ID_vec <- sub("^.*([0-9]{6})$", "\\1", raw_id)
    BV_vec <- suppressWarnings(as.numeric(dfB[[2]][101:600]))

    params <- tuned_params_by_ns[[as.character(ns)]]

    if (is.null(params)) {
      stop(sprintf("No tuned parameters were defined for ns=%d", ns))
    }

    # Push data into the AMPL model and solve with Gurobi through AMPL.
    optorch_set_data(ampl, b = BV_vec, c = C_mat, l = L, u = U, ns = ns)

    res <- optorch_solve(
      ampl,
      gurobi_params = params,
      log_file = file.path(iter_dir, sprintf("gurobi_rep_%02d_ns_%02d.log", iter, ns)),
      time_limit = time_limit
    )

    status <- res$status
    objval <- res$objval
    runtime <- res$runtime

    if (identical(status, "OPTIMAL")) {

      p_vals <- res$p
      r_vals <- res$r

      selected_idx <- which(r_vals > 0.5)
      n_selected <- length(selected_idx)

      cat(sprintf(
        "ns=%d iter=%d status=OK obj=%.4f n=%d runtime=%.2f\n",
        ns, iter, objval, n_selected, runtime
      ))

      saveRDS(
        p_vals,
        file.path(iter_dir, sprintf("p_rep_%02d_ns_%02d.rds", iter, ns))
      )

      saveRDS(
        r_vals,
        file.path(iter_dir, sprintf("r_rep_%02d_ns_%02d.rds", iter, ns))
      )

      row_id_iter <- row_id_iter + 1
      results_summary_iter[[row_id_iter]] <- data.frame(
        iter = iter,
        ns = ns,
        status = status,
        objval = objval,
        runtime = runtime,
        n_selected = n_selected,
        solve_result = res$solve_result,
        solve_message = res$solve_message,
        stringsAsFactors = FALSE
      )

      if (n_selected > 0) {

        sel_ids <- ID_vec[selected_idx]
        sel_bv <- BV_vec[selected_idx]

        sel_df <- data.frame(
          ID = sel_ids,
          iter = iter,
          ns = ns,
          bv = sel_bv,
          p = p_vals[selected_idx],
          r = r_vals[selected_idx],
          stringsAsFactors = FALSE
        )

        sel_row_id_iter <- sel_row_id_iter + 1
        selected_rows_iter[[sel_row_id_iter]] <- sel_df

        C_sub <- C_mat[selected_idx, selected_idx, drop = FALSE]

        write.csv(
          C_sub,
          file.path(iter_dir, sprintf("Csub_rep_%02d_ns_%02d.csv", iter, ns)),
          row.names = TRUE
        )
      }

    } else {

      cat(sprintf("ns=%d iter=%d status=%s (%s)\n", ns, iter, status, res$solve_message))

      row_id_iter <- row_id_iter + 1
      results_summary_iter[[row_id_iter]] <- data.frame(
        iter = iter,
        ns = ns,
        status = status,
        objval = objval,
        runtime = runtime,
        n_selected = NA_integer_,
        solve_result = res$solve_result,
        solve_message = res$solve_message,
        stringsAsFactors = FALSE
      )
    }

    if (length(results_summary_iter) > 0) {
      write.csv(
        do.call(rbind, results_summary_iter),
        file.path(iter_dir, "summary.csv"),
        row.names = FALSE
      )
    }

    if (length(selected_rows_iter) > 0) {
      write.csv(
        do.call(rbind, selected_rows_iter),
        file.path(iter_dir, "selected_individuals.csv"),
        row.names = FALSE
      )
    }
  }
}

cat("Done. Results are grouped by ns, then per iteration.\n")
