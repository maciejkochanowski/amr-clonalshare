# The R packages of the comparison with published procedures, at the versions
# of its first run, installed into the library named by R_LIBS_BENCH:
#
#   R_LIBS_BENCH=$ROOT/rlib Rscript benchmarks/competitors/install_r_bench.R
#
# The build of R 4.3.3 on the cluster of the revision campaign lacks the
# recommended packages, so they are installed first at the versions R 4.3.3
# ships with.
lib <- Sys.getenv("R_LIBS_BENCH")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
for (p in c("remotes", "BiocManager"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib)
have <- function(p) tryCatch(format(packageVersion(p, lib.loc = .libPaths())), error = function(e) "")
want <- c(MASS = "7.3-60.0.1", lattice = "0.22-5", Matrix = "1.6-5", nlme = "3.1-164", mgcv = "1.9-1",
          survival = "3.5-8", interval = "1.1-1.0", rankICC = "1.0.2", lme4 = "1.1-35.1", coin = "1.4-3")
for (p in names(want)) {
  if (p == "interval" && have("Icens") == "")
    BiocManager::install("Icens", lib = lib, update = FALSE, ask = FALSE)
  if (have(p) != gsub("-", ".", want[[p]]))
    remotes::install_version(p, version = want[[p]], lib = lib, upgrade = "never", dependencies = NA)
}
cat(R.version.string, "\n")
for (p in names(want)) cat(p, format(packageVersion(p)), "\n")
