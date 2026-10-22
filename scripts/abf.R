suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
config <- fromJSON(args[1])
stopifnot(as.character(packageVersion("coloc")) == config$software$coloc)
task <- fromJSON(args[2])

dataset <- function(d) {
  stopifnot(isTRUE(d$scale_certified), all(is.finite(d$beta)),
            all(is.finite(d$se)), all(d$se > 0))
  x <- list(beta = d$beta, varbeta = d$se^2, snp = task$snp,
            position = task$position, type = d$type)
  if (d$type == "cc") stopifnot(d$effect_scale == "log_odds_ratio")
  if (d$type == "quant") {
    x$N <- d$N
    x$MAF <- d$MAF
  }
  x
}

fit <- coloc.abf(dataset(task$gwas), dataset(task$qtl),
                 p1 = config$abf$p1, p2 = config$abf$p2,
                 p12 = config$abf$p12)
stopifnot(fit$summary[["nsnps"]] == length(task$snp))
write_json(list(task_id = task$task_id, n_variants = length(task$snp),
                pp_h4 = unname(fit$summary[["PP.H4.abf"]])),
           args[3], auto_unbox = TRUE, pretty = TRUE, digits = 17)
