# optorch_ampl.R
#
# R interface to the OptOrch AMPL models through rAMPL.
#
# The AMPL model files in 01_models/ are the single source of truth for the
# optimization formulation. This file only moves data between R and AMPL,
# configures the solver, and collects results. It replaces the direct
# gurobi-R-API model construction (sparse A matrix, quadcon list, ...) that the
# original reproduction scripts and the Shiny app used.
#
# Requirements:
#   - rAMPL (https://rampl.ampl.com)
#   - An AMPL installation. Either put the folder containing ampl(.exe) on the
#     system PATH or set the environment variable AMPL_HOME to that folder.
#   - A solver that handles non-convex MIQCP. The manuscript used Gurobi, which
#     is the default here. Any AMPL solver name can be passed as `solver`.
#
# Typical use:
#   source("01_models/optorch_ampl.R")
#   ampl <- optorch_ampl_start()
#   optorch_load_model(ampl, model_dir = "01_models", scenario = "status_number")
#   optorch_set_data(ampl, b = BV_vec, c = C_mat, l = 0.01, u = 0.15, ns = 30)
#   res <- optorch_solve(ampl, gurobi_params = list(MIPGap = 0.005, NonConvex = 2))
#   res$status; res$objval; res$p; res$r
#   ampl$close()

library(rAMPL)

#---------------------------------------------------------------
# AMPL session ---------------------------------------------------

#' Start an AMPL session.
#'
#' @param ampl_dir Folder containing the AMPL executable. Defaults to the
#'   AMPL_HOME environment variable; if that is empty, AMPL is looked up on PATH.
#' @param solver AMPL solver name (default "gurobi").
#' @return An rAMPL AMPL object.
optorch_ampl_start <- function(ampl_dir = Sys.getenv("AMPL_HOME", ""), solver = "gurobi") {
  ampl <- if (nzchar(ampl_dir)) {
    new(AMPL, new(Environment, normalizePath(ampl_dir, winslash = "/", mustWork = TRUE)))
  } else {
    new(AMPL)
  }
  ampl$setOption("solver", solver)
  ampl
}

#---------------------------------------------------------------
# Model ----------------------------------------------------------

#' Load the OptOrch AMPL formulation for a scenario.
#'
#' scenario = "status_number"        : OptOrch_core.mod (manuscript Scenario 1)
#' scenario = "pollen_contamination" : OptOrch_core.mod + OptOrch_supp6.mod
#'                                     (manuscript Scenario 2). supp6 replaces
#'                                     the total-contribution (c2) and
#'                                     status-number (c1) constraints of the
#'                                     core model, so those two are dropped.
#'
#' @param ampl An AMPL object from optorch_ampl_start().
#' @param model_dir Path to the 01_models folder.
#' @param scenario One of "status_number", "pollen_contamination".
optorch_load_model <- function(ampl, model_dir, scenario = c("status_number", "pollen_contamination")) {
  scenario <- match.arg(scenario)
  model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)

  ampl$reset()
  ampl$read(file.path(model_dir, "001_core", "OptOrch_core.mod"))

  if (scenario == "pollen_contamination") {
    ampl$read(file.path(model_dir, "002_tailored", "OptOrch_supp6.mod"))
    ampl$eval("drop c1; drop c2;")
  }

  invisible(scenario)
}

#---------------------------------------------------------------
# Data -----------------------------------------------------------

#' Push candidate data into the loaded AMPL model.
#'
#' @param ampl An AMPL object with a model loaded.
#' @param b Numeric vector of breeding values (length nc).
#' @param c Numeric nc x nc coancestry matrix (already on the coancestry scale,
#'   i.e. G/2, and already made PSD if desired).
#' @param l,u Lower/upper contribution bounds: a scalar or a length-nc vector.
#' @param ns Target status number.
#' @param contamination_rate,external_pollen_parents Only used by the
#'   pollen-contamination scenario (parameters C and nd of OptOrch_supp6.mod).
optorch_set_data <- function(ampl, b, c, l, u, ns,
                             contamination_rate = NULL,
                             external_pollen_parents = NULL) {
  nc <- length(b)
  if (!all(dim(c) == c(nc, nc))) {
    stop(sprintf("c must be a %d x %d matrix.", nc, nc), call. = FALSE)
  }
  if (length(l) == 1) l <- rep(l, nc)
  if (length(u) == 1) u <- rep(u, nc)
  if (length(l) != nc || length(u) != nc) {
    stop("l and u must be scalars or vectors of length nc.", call. = FALSE)
  }

  ampl$getParameter("nc")$set(nc)
  ampl$getParameter("ns")$set(ns)

  # 1-index parameters over S = 1..nc: values are assigned in set order.
  ampl$getParameter("b")$setValues(as.numeric(b))
  ampl$getParameter("l")$setValues(as.numeric(l))
  ampl$getParameter("u")$setValues(as.numeric(u))

  # 2-index parameter c{S,S}: long data.frame (i, j, value), column-major.
  c_df <- data.frame(
    i = rep.int(seq_len(nc), times = nc),
    j = rep(seq_len(nc), each = nc),
    c = as.numeric(c)
  )
  ampl$getParameter("c")$setValues(c_df)

  if (!is.null(contamination_rate)) {
    ampl$getParameter("C")$set(contamination_rate)
  }
  if (!is.null(external_pollen_parents)) {
    ampl$getParameter("nd")$set(external_pollen_parents)
  }

  invisible(nc)
}

#---------------------------------------------------------------
# Solver options -------------------------------------------------

#' Build the gurobi_options string for the AMPL Gurobi driver.
#'
#' Gurobi parameters are passed with their native Gurobi names through the
#' driver's `optionnative` option, so the tuned parameter sets from the
#' manuscript can be reused verbatim (e.g. list(MIPGap = 0.005, NonConvex = 2)).
#'
#' @param gurobi_params Named list of native Gurobi parameters.
#' @param outlev 1 to print the Gurobi log to the console, 0 to stay silent.
#' @param log_file Optional path for the Gurobi log file.
#' @param time_limit Optional time limit in seconds.
optorch_gurobi_options <- function(gurobi_params = list(), outlev = 1,
                                   log_file = NULL, time_limit = NULL) {
  # OutputFlag/LogFile/TimeLimit are driver-level options in AMPL; translate
  # them if they show up in a legacy parameter list.
  if (!is.null(gurobi_params$OutputFlag)) {
    outlev <- as.integer(gurobi_params$OutputFlag)
    gurobi_params$OutputFlag <- NULL
  }
  if (!is.null(gurobi_params$LogFile)) {
    log_file <- gurobi_params$LogFile
    gurobi_params$LogFile <- NULL
  }
  if (!is.null(gurobi_params$TimeLimit)) {
    time_limit <- gurobi_params$TimeLimit
    gurobi_params$TimeLimit <- NULL
  }

  opts <- sprintf("outlev=%d", as.integer(outlev))
  if (!is.null(time_limit)) {
    opts <- c(opts, sprintf("timelim=%s", format(time_limit, scientific = FALSE)))
  }
  if (!is.null(log_file)) {
    opts <- c(opts, sprintf("logfile='%s'", normalizePath(log_file, winslash = "/", mustWork = FALSE)))
  }
  for (nm in names(gurobi_params)) {
    opts <- c(opts, sprintf('optionnative="%s %s"', nm,
                            format(gurobi_params[[nm]], scientific = FALSE)))
  }
  paste(opts, collapse = " ")
}

#---------------------------------------------------------------
# Solve ----------------------------------------------------------

#' Solve the loaded OptOrch model and collect results.
#'
#' @inheritParams optorch_gurobi_options
#' @return A list with
#'   status        "OPTIMAL" when AMPL reports solve_result == "solved";
#'                 otherwise the upper-cased AMPL solve_result ("LIMIT",
#'                 "INFEASIBLE", ...). Kept compatible with the Gurobi status
#'                 strings used in the original output files.
#'   solve_result  AMPL solve_result string.
#'   solve_message Solver message returned to AMPL.
#'   objval        Objective value (NA unless a solution is available).
#'   runtime       Solver wall-clock time in seconds.
#'   p, r          Numeric vectors (length nc) in candidate order, or NULL.
optorch_solve <- function(ampl, gurobi_params = list(), outlev = 1,
                          log_file = NULL, time_limit = NULL) {
  ampl$setOption("gurobi_options",
                 optorch_gurobi_options(gurobi_params, outlev, log_file, time_limit))

  t0 <- proc.time()[["elapsed"]]
  ampl$solve()
  elapsed_r <- proc.time()[["elapsed"]] - t0

  solve_result <- as.character(ampl$getValue("solve_result"))
  # Keep only the first line of the solver message (e.g.
  # "Gurobi 13.0.2: optimal solution; objective 1.51"); the remaining lines
  # are iteration counts that are already in the log file.
  solve_message <- trimws(strsplit(as.character(ampl$getValue("solve_message")), "\n")[[1]][1])
  runtime <- suppressWarnings(as.numeric(ampl$getValue("_solve_elapsed_time")))
  if (is.na(runtime)) runtime <- elapsed_r

  status <- if (identical(solve_result, "solved")) "OPTIMAL" else toupper(solve_result)

  # A solution vector is available when solved, or when a limit was hit with
  # a feasible incumbent (AMPL reports "limit" and keeps the incumbent).
  has_solution <- solve_result %in% c("solved", "limit") &&
    !grepl("no (feasible )?solution|infeasible", solve_message, ignore.case = TRUE)

  objval <- NA_real_
  p <- NULL
  r <- NULL
  if (has_solution) {
    objval <- ampl$getObjective("gain")$value()
    p_df <- ampl$getVariable("p")$getValues()
    r_df <- ampl$getVariable("r")$getValues()
    p <- p_df[[2]][order(p_df[[1]])]
    r <- r_df[[2]][order(r_df[[1]])]
  }

  list(
    status = status,
    solve_result = solve_result,
    solve_message = solve_message,
    objval = objval,
    runtime = runtime,
    p = p,
    r = r
  )
}

#---------------------------------------------------------------
# Convenience: one-shot run --------------------------------------

#' Load model, set data, solve. Convenience wrapper for single runs.
optorch_run <- function(ampl, model_dir, scenario, b, c, l, u, ns,
                        contamination_rate = NULL, external_pollen_parents = NULL,
                        gurobi_params = list(), outlev = 1, log_file = NULL,
                        time_limit = NULL) {
  optorch_load_model(ampl, model_dir, scenario)
  optorch_set_data(ampl, b = b, c = c, l = l, u = u, ns = ns,
                   contamination_rate = contamination_rate,
                   external_pollen_parents = external_pollen_parents)
  optorch_solve(ampl, gurobi_params = gurobi_params, outlev = outlev,
                log_file = log_file, time_limit = time_limit)
}
