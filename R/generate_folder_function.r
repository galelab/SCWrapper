#' generate folder function
#'
#' This generates new folder (internal function only)
#' @param foldername folder to genrate
#' @keywords remove samples
#' @export

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}