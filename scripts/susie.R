suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
config <- fromJSON(args[1])
stopifnot(as.character(packageVersion("coloc")) == config$software$coloc,
          as.character(packageVersion("susieR")) == config$software$susieR,
          config$susie$maxit == 2000)
task <- fromJSON(args[2])
n <- length(task$snp)

dataset <- function(d) {
  stopifnot(d$N > 1, isTRUE(d$scale_certified), identical(d$snp, task$snp))
  con <- file(args[3], "rb")
  values <- readBin(con, "double", n = n * n, size = 8, endian = "little")
  close(con)
  stopifnot(length(values) == n * n, all(is.finite(values)))
  ld <- matrix(values, n, n, byrow = TRUE,
               dimnames = list(task$snp, task$snp))
  stopifnot(max(abs(ld - t(ld))) < 1e-8,
            max(abs(diag(ld) - 1)) < 1e-8,
            min(eigen(ld, symmetric = TRUE, only.values = TRUE)$values) >= -1e-8)
  x <- list(beta = d$beta, varbeta = d$se^2, snp = task$snp,
            position = task$position, N = d$N, LD = ld, type = d$type)
  if (d$type == "quant") {
    if (!is.null(d$sdY)) x$sdY <- d$sdY else x$MAF <- d$MAF
  }
  x
}

gwas <- dataset(task$gwas)
qtl <- dataset(task$qtl)
fit <- function(d) {
  set.seed(config$susie$seed)
  coloc::runsusie(d, L = config$susie$L, maxit = config$susie$maxit,
                  tol = config$susie$tol, coverage = config$susie$coverage,
                  min_abs_corr = config$susie$min_abs_corr,
                  repeat_until_convergence = config$susie$repeat_until_convergence,
                  estimate_residual_variance = config$susie$estimate_residual_variance,
                  z_ld_weight = config$susie$z_ld_weight)
}
run_once <- function() {
  a <- fit(gwas)
  b <- fit(qtl)
  stopifnot(isTRUE(a$converged), isTRUE(b$converged),
            all(is.finite(a$pip)), all(is.finite(b$pip)))
  if (!length(a$sets$cs) || !length(b$sets$cs))
    return(list(status = "NOT_EVALUABLE", reason = "NO_VALID_SIGNAL_PAIR"))
  result <- coloc::coloc.susie(a, b, p1 = config$abf$p1,
                                p2 = config$abf$p2, p12 = config$abf$p12)
  pairs <- as.data.frame(result$summary)
  posterior <- as.matrix(pairs[, paste0("PP.H", 0:4, ".abf"), drop = FALSE])
  stopifnot(all(is.finite(posterior)), all(posterior >= 0 & posterior <= 1),
            all(abs(rowSums(posterior) - 1) < 1e-8))
  list(status = "SUCCESS", pp_h4 = max(pairs$PP.H4.abf),
       pair_summary = pairs, gwas_iterations = a$niter,
       qtl_iterations = b$niter, gwas_pip = a$pip, qtl_pip = b$pip)
}

first <- run_once()
second <- run_once()
stopifnot(isTRUE(all.equal(first, second, tolerance = 1e-12)))
first$task_id <- task$task_id
first$max_iter <- config$susie$maxit
first$repeatability <- "PASS_IDENTICAL_SEEDED_REPEAT"
write_json(first, args[4], auto_unbox = TRUE, pretty = TRUE, digits = 17)
