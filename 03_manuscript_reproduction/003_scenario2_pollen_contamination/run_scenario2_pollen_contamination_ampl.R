library(rAMPL)
library(data.table)

# rAMPL version of run_scenario2_pollen_contamination.R.
#
# The optimization model is no longer assembled in R for the Gurobi API.
# Instead, the AMPL formulation (01_models/001_core/OptOrch_core.mod with the
# pollen-contamination constraints of 01_models/002_tailored/OptOrch_supp6.mod)
# is loaded through rAMPL (see 01_models/optorch_ampl.R), the data are pushed
# from R, and Gurobi is called through AMPL with the same tuned parameters.
# Output files have the same names and columns as the original script, plus
# two extra columns in summary.csv (solve_result, solve_message).

#---------------------------------------------------------------
# 0) Path ------------------------------------------------------

# Set the working directory to:
# 03_manuscript_reproduction/003_scenario2_pollen_contamination/

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
#   OPTORCH_NS_GRID=30 OPTORCH_ITER_GRID=1 OPTORCH_TIME_LIMIT=120 Rscript run_scenario2_pollen_contamination.R
if (nzchar(Sys.getenv("OPTORCH_NS_GRID"))) {
  ns_grid <- as.integer(strsplit(Sys.getenv("OPTORCH_NS_GRID"), ",")[[1]])
}
if (nzchar(Sys.getenv("OPTORCH_ITER_GRID"))) {
  iter_grid <- as.integer(strsplit(Sys.getenv("OPTORCH_ITER_GRID"), ",")[[1]])
}
time_limit <- if (nzchar(Sys.getenv("OPTORCH_TIME_LIMIT"))) as.numeric(Sys.getenv("OPTORCH_TIME_LIMIT")) else NULL

# Pollen contamination
nd <- 100
C <- 0.3

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

# Core model + OptOrch_supp6.mod (effective contribution 1 - C/2 and the
# contamination-adjusted status-number constraint replace c1 and c2).
optorch_load_model(ampl, MODEL_DIR, scenario = "pollen_contamination")

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

    ext_path <- file.path(DATA_DIR, sprintf("rep%d_Np2_Data.csv", iter))

    if (is.na(c_path) || is.na(b_path) || !file.exists(ext_path)) {
      cat(sprintf(
        "ns=%d iter=%d skipped because one or more input files are missing\n",
        ns, iter
      ))
      next
    }

    dfExt <- fread(ext_path)

    if (!("BV" %in% names(dfExt))) {
      stop(sprintf("Column 'BV' was not found in external pollen data: %s", ext_path))
    }

    BV_ext_mean <- mean(as.numeric(dfExt$BV), na.rm = TRUE)

    # C matrix: rows and columns 102-601 in the original CSV
    dfC <- read.csv(
      c_path,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    C_block <- dfC[102:601, 102:601, drop = FALSE]
    C_mat <- apply(C_block, 2, function(x) suppressWarnings(as.numeric(x)))

    # Convert the genomic relationship matrix G to the coancestry matrix C.
    C_mat <- as.matrix(C_mat) * 0.5

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

    # Push data into the AMPL model (including C and nd for OptOrch_supp6.mod)
    # and solve with Gurobi through AMPL.
    optorch_set_data(
      ampl,
      b = BV_vec, c = C_mat, l = L, u = U, ns = ns,
      contamination_rate = C,
      external_pollen_parents = nd
    )

    res <- optorch_solve(
      ampl,
      gurobi_params = params,
      log_file = file.path(iter_dir, sprintf("gurobi_rep_%02d_ns_%02d.log", iter, ns)),
      time_limit = time_limit
    )

    status <- res$status
    objval <- res$objval
    runtime <- res$runtime

    if (!identical(status, "OPTIMAL")) {

      write.csv(
        data.frame(
          iter = iter,
          ns = ns,
          status = status,
          objval = objval,
          runtime = runtime,
          solve_result = res$solve_result,
          solve_message = res$solve_message,
          stringsAsFactors = FALSE
        ),
        file.path(iter_dir, "summary.csv"),
        row.names = FALSE
      )

      cat(sprintf("NOT OPTIMAL ns=%d iter=%d status=%s (%s)\n", ns, iter, status, res$solve_message))
      next
    }

    p_vals <- res$p
    r_vals <- res$r

    Gain_internal <- sum(p_vals * BV_vec)
    Gain_total <- Gain_internal + (C / 2) * BV_ext_mean

    cat(sprintf(
      "ns=%d iter=%d status=OK obj=%.4f runtime=%.2f sum(p)=%.6f Gain_internal=%.6f BVext=%.6f Gain_total=%.6f\n",
      ns, iter, objval, runtime, sum(p_vals), Gain_internal, BV_ext_mean, Gain_total
    ))

    saveRDS(
      p_vals,
      file.path(iter_dir, sprintf("p_rep_%02d_ns_%02d.rds", iter, ns))
    )

    saveRDS(
      r_vals,
      file.path(iter_dir, sprintf("r_rep_%02d_ns_%02d.rds", iter, ns))
    )

    write.csv(
      data.frame(
        iter = iter,
        ns = ns,
        status = status,
        objval = objval,
        runtime = runtime,
        n_selected = sum(r_vals > 0.5),
        sum_p = sum(p_vals),
        BV_ext_mean = BV_ext_mean,
        Gain_internal = Gain_internal,
        Gain_total = Gain_total,
        solve_result = res$solve_result,
        solve_message = res$solve_message,
        stringsAsFactors = FALSE
      ),
      file.path(iter_dir, "summary.csv"),
      row.names = FALSE
    )

    selected_idx <- which(r_vals > 0.5)

    if (length(selected_idx) > 0) {

      write.csv(
        data.frame(
          ID = ID_vec[selected_idx],
          iter = iter,
          ns = ns,
          bv = BV_vec[selected_idx],
          p = p_vals[selected_idx],
          r = r_vals[selected_idx],
          stringsAsFactors = FALSE
        ),
        file.path(iter_dir, "selected_individuals.csv"),
        row.names = FALSE
      )

      C_sub <- C_mat[selected_idx, selected_idx, drop = FALSE]

      write.csv(
        C_sub,
        file.path(iter_dir, sprintf("Csub_rep_%02d_ns_%02d.csv", iter, ns)),
        row.names = TRUE
      )
    }
  }
}

cat("Done.\n")
