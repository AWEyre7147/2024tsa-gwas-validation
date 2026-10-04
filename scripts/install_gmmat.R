args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1])
folder <- file.path(root, "data", "tools", "association")
dir.create(folder, recursive=TRUE, showWarnings=FALSE)
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=300)
install_pinned <- function(package, version) {
    archive <- file.path(folder, paste0(package,"_",version,".tar.gz"))
    if (!file.exists(archive)) {
        current <- paste0("https://cloud.r-project.org/src/contrib/",basename(archive))
        historical <- paste0("https://cloud.r-project.org/src/contrib/Archive/",package,"/",basename(archive))
        tryCatch(download.file(current, archive, mode="wb"),
            error=function(e) download.file(historical, archive, mode="wb"))
    }
    if (!requireNamespace(package, quietly=TRUE) || packageVersion(package) != version)
        install.packages(archive, repos=NULL, type="source")
    stopifnot(requireNamespace(package, quietly=TRUE), packageVersion(package)==version)
}
install_pinned("CompQuadForm", "1.4.4")
stopifnot(requireNamespace("data.table",quietly=TRUE),packageVersion("data.table")=="1.18.6.1")
archive <- file.path(folder, "GMMAT_1.5.0.tar.gz")
if (!file.exists(archive)) {
    tryCatch(download.file("https://cloud.r-project.org/src/contrib/GMMAT_1.5.0.tar.gz", archive, mode="wb"),
        error=function(e) download.file("https://cloud.r-project.org/src/contrib/Archive/GMMAT/GMMAT_1.5.0.tar.gz",archive,mode="wb"))
}
if (!requireNamespace("GMMAT", quietly=TRUE) || packageVersion("GMMAT") != "1.5.0") {
    install.packages(archive, repos=NULL, type="source")
}
stopifnot(requireNamespace("GMMAT", quietly=TRUE), packageVersion("GMMAT") == "1.5.0")
writeLines(capture.output(sessionInfo()), file.path(folder, "R-session.private.txt"))
cat("GMMAT", as.character(packageVersion("GMMAT")), "installed successfully.\n")
