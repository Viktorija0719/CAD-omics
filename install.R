# CRAN packages
install.packages(c(
  "tidyr",
  "readr",
  "caret"
), repos = "https://cran.r-project.org")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cran.r-project.org")

BiocManager::install(c(
  "edgeR",
  "limma"
), ask = FALSE, update = TRUE)

# Optional: validate Bioconductor packages
BiocManager::valid()
