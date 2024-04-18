
# Utility functions

#' Append console output if verbose
#'
#' @param ... function to pass to ca
#' @return console output to outfile
#' @keywords internal
vcat <- function(..., append = TRUE) if (verbose) cat(..., file = outfile, append = append)

#' Add carriage return to .csv files
#'
#' @param dataDir directory containing .csv files
#' @description Adds a carriage return (end of line character) to .csv files. Note that this is recursive by default. This behavior can be changed by setting "recursive = FALSE".
#' @return modified .csv files
#' @keywords internal
crAdd <- function(dataDir, ..) {
  ls <- list.files(path = dataDir, pattern = ".csv$", recursive = TRUE, full.names = TRUE)
  sapply(ls, function(file) write.table(
    "", file = file, sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE))
}


#' Check whether an object contains SWEDDIE core tables
#'
#' @param x The object to check
#' @description A valid SWEDDIE object is a list with the following elements,
#' all of which must be \code{\link{data.frame}} objects:
#' \itemize{
#' \item{\code{experiment}}{ Experiment}
#' \item{\code{site}}{ Site}
#' \item{\code{plot}}{ Plot}
#' }
#' @return TRUE or FALSE.
#' @keywords internal
check_sweddie_core <- function(x) {
  # Database is a list and must have all the following data frames
  tables <- c(
    "experiment", "site", "plot"
  )
  is.list(x) &&
    identical(sort(tables), sort(names(x))) &&
    all(sapply(x, class) == "data.frame")
}


#' Check that column names in the template and meta data files match.
#'
#' @param metadata list of core template files
#' @param tableName specific table name to check
#' @param datIn input data to check
#' @param err error counter; defaults to 0
#' @return error counter "err"
#' @note This is typically called only from \code{\link{checkTemplateFiles}}.
#' @keywords internal
# checks column names
checkColNms <- function(metadata, tableName, datIn, err = 0, ...) {

  # check for template
  if (missing(metadata)) metadata <- readMeta()

  # check tableNames
  if (!names(datIn)[which(names(datIn) == tableName)] %in% names(metadata)[which(names(metadata) == tableName)])  {
    err <- err + 1
    vcat(warning(
      paste("\t", names(datIn)[tableName], "is missing in data file\n")))
  }

  # check for missing cols in data
  miss <- setdiff(colnames(metadata[[tableName]]), colnames(datIn[[tableName]]))

  # check for extra cols in data
  xtra <- setdiff(colnames(datIn[[tableName]]), colnames(metadata[[tableName]]))

  if (length(miss) > 0 | length(xtra) > 0) {
    err <- err + 1
    vcat(warning(
      paste("\t", names(datIn[tableName]), "table names do not match template\n")))
    if (length(miss) > 0) vcat("\t\tColumn names missing:", miss, "\n")
    if (length(xtra > 0)) vcat("\t\tColumn names extra:", xtra, "\n")
  }
  return(err)
}

#' Check that a column is strictly numeric.
#'
#' @param x Column values, a vector
#' @param xname Column name
#' @return Nothing (run for its warning side effect).
#' @keywords internal
check_numeric <- function(x, xname) {
  stopifnot(is.character(xname))

  if (!is.numeric(type.convert(x, as.is = FALSE))) {
    warning("Non-numeric values in ", xname, " column")
  }
}

#' Read template files
#'
#' @param metaDir core template directory
#' @param return should an empty template list be returned ("template") or the meta data file itself?
#' @return A list with core template tables or core meta data
#' @keywords internal
readMeta <- function(metaDir, return = "template") {
  # get dir name as needed
  if (missing(metaDir)) {
    metaDir <- "../data/sweddie/metadata"
  }
  # create template list
  ls <- list.files(metaDir, full.names = TRUE)
  ls <- ls[grep(paste(c("exp", "sit", "plt"), collapse = "|"), ls)]
  meta <- lapply(
    setNames(ls, nm = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(ls))),
    read.csv)
  if (return == "template") {
    lapply(
      setNames(meta, nm = c("experiment", "plot", "site")), function(x) {
        setNames(data.frame(matrix(ncol = nrow(x), nrow = 0)), x[ , 1])
      })
  } else {
    return(meta)
  }
}
