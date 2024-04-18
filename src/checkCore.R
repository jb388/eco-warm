#' QAQC
#'
#' @description Checks template files for data coherence, formatting, and data entry errors
#'
#' @details This function can also be called from the \href{https://soilradiocarbon.org}{ISRaD website}.
#' @param file File path for template file to be checked
#' @param writeQCreport If TRUE, a text report of the QC output will be written to the outfile. Default is FALSE
#' @param outfile_QAQC Filename of the output file (if writeQCreport is TRUE). Default is NULL, with the outfile being written to the directory where the template file is stored and named according to the file being checked.
#' @param summaryStats Prints summary statistics. Default is TRUE.
#' @param dataReport Prints list structure of database. Default is FALSE.
#' @param checkdoi Set to FALSE if you do not want the QAQC check to validate DOIs (if TRUE this will be time consuming). Default is TRUE.
#' @param verbose Set to TRUE to print results of function to console. Default is TRUE.
#' @param local Set to FALSE to fetch most up-to-date template and template info files. If TRUE, the local files or files from CRAN package will be used. Default is TRUE.
#' @import dplyr
#' @importFrom RCurl url.exists
#' @importFrom readxl read_excel excel_sheets
#' @importFrom httr HEAD
#' @importFrom rio import
#' @importFrom utils type.convert
#' @export
#' @examples
#' \donttest{
#' # Load example dataset Gaudinski_2001
#' entry <- ISRaD::Gaudinski_2001
#' # Save as .xlsx file
#' ISRaD.save.entry(
#'   entry = entry,
#'   template_file = system.file("extdata", "ISRaD_Master_Template.xlsx", package = "ISRaD"),
#'   outfile = file.path(tempdir(), "Gaudinski_2001.xlsx")
#' )
#' # Run QAQC
#' QAQC(file.path(tempdir(), "Gaudinski_2001.xlsx"))
#' }
#'
checkCore <- function(file, writeQCreport = FALSE, outfile_QAQC = "", summaryStats = TRUE,
                 dataReport = FALSE, checkdoi = TRUE, verbose = TRUE, local = TRUE) {
  stopifnot(is.character(file))
  stopifnot(is.logical(writeQCreport))
  stopifnot(is.character(outfile_QAQC))
  stopifnot(is.logical(summaryStats))
  stopifnot(is.logical(dataReport))
  stopifnot(is.logical(checkdoi))
  stopifnot(is.logical(verbose))

  vcat <- function(..., append = TRUE) if (verbose) cat(..., file = outfile_QAQC, append = append)
  
  # start error count at 0
  error <- 0
  # start note count at 0
  note <- 0

  if (writeQCreport) {
    if (outfile_QAQC == "") {
      outfile_QAQC <- file.path(dirname(file), "QAQC", paste0("QAQC_", gsub("\\.xlsx", ".txt", basename(file))))
    }
  }

  vcat("\nFile:", basename(file))
  # message("\nTime:", as.character(Sys.time()), "\n", file=outfile_QAQC, append = TRUE)

  ##### check file extension #####
  vcat("\n\nChecking file type...")
  if (!grep(".xlsx", file) == 1) {
    vcat("\nWARNING: ", file, " is not the current file type (should have '.xlsx' extension)")
    error <- error + 1
  }


