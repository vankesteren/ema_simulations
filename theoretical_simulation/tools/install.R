# install script to make sure the correct versions of the packages are installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# regsem for xmed
devtools::install_cran("Rsolnp")
if (!requireNamespace("regsem") || packageVersion("regsem") != "1.0.6")
  devtools::install_version("regsem", "1.0.6")

# pbapply for parallel processing
devtools::install_cran("pbapply", quiet = TRUE)

# HIMA
devtools::install_github(
  "YinanZheng/HIMA", ref = "d9724d3c78502d66cdf9cf7130c74e8ebb9cefaf"
)

# CMF
devtools::install_github(
  "vankesteren/cmfilter", ref = "0f016fd8d3a25b068e9fa97ad9f658c61cf186ab"
)
