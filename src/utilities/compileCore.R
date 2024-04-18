#' Compile SWEDDIE core files
#'
#' Compiles experiment, site, and plot data
#'
#' @param DIR Directory where completed QAQCed template files are stored. Defaults to "./data/swaedie"
#' @param write_report Boolean flag to write a log file of the
#' compilation. File will be in the specified
#' dataset_directory at "database/coreLog_yyyymmddtttt.txt" (string after underscore is an explicit time stamp).
#' @param write_out Set to TRUE to write compiled core data to output file. If FALSE, function will simply return the database as a list.
#' @param output_ext Set file type for output file.
#' @param verbose Set to TRUE to print results of function to console.
#'
#' @export
#' @importFrom dplyr mutate_all setdiff
#' @examples
compileCore <- function(
    DIR = "../data/sweddie",
    write_report = TRUE,
    write_out = TRUE,
    output_ext = c(".csv", ".rda"),
    verbose = TRUE) {

  # Check inputs
  stopifnot(dir.exists(DIR))
  stopifnot(is.logical(write_out))
  stopifnot(is.logical(write_report))
  stopifnot(is.character(output_ext))
  stopifnot(is.logical(verbose))

  # Constants
  DB_DIR <- "database"
  S_DIR <- "siteData"
  LIST_FILE <- "coreData.rda"
  TIMESTAMP <- format(Sys.time(), "%y%m%d-%H%M")

  # Set output file
  outfile <- ""
  if (write_report) {
    outfile <- file.path(DIR, DB_DIR, paste0("logs/coreLog", "_", TIMESTAMP, ".txt"))
    file.create(outfile)
  }

  # Start writing in the output file
  vcat("SWAEDIE Compilation Log \n",
       "\n", as.character(Sys.time()),
       "\n", rep("-", 15), "\n")


  # Get the tables stored in the template sheets
  metadata <- readMeta()

  vcat("\n\nCompiling data files in", S_DIR, "\n", rep("-", 30), "\n")

  data_dirs <- list.dirs(file.path(DIR, S_DIR), full.names = TRUE, recursive = FALSE)
  if (!length(data_dirs)) {
    vcat(warning("No data directories found!\n"))
    return(NULL)
  }

  # ensure EOL carriage return present
  invisible(crAdd(file.path(DIR, S_DIR)))

  vcat("Compiling and checking site data...\n")
    #   pb <- txtProgressBar(min = 0, max = length(data_dirs), style = 3)
    # }
    #
    # check if previous database object exists in database directory, and only update file if new data exisit
  if (file.exists(file.path(DIR, DB_DIR, LIST_FILE))) {

    # load existing database
    load(file.path(dataset_directory, DB_DIR, LIST_FILE)) # obj "coreDat"

    # convert to character and coerce to list of data frames
    coreDat_chr <- lapplydf(lapply(coreData, function(x) lapply(x, as.character)))

    # remove old version of ISRaD
    rm(coreDat)

    # Split each table by entry_name
    coreDat_old <- lapply(coreDat_chr, function(x) split(x, x$exp_name))
  } else {
    database <- metadata
  }

  # compile new templates and check against existing data
  for (d in seq_along(data_dirs)) {

    # get site name
    site_nm <- basename(data_dirs[d])

    vcat(paste(site_nm, "\n"))

    # get tables
    tbls <- list.files(data_dirs[d], full.names = TRUE)

    # check that required files are present
    ix <- which(!names(metadata) %in% sub("\\..*", "", basename(tbls)))
    if (length(ix) > 0) {
      vcat(warning(
        paste0(paste(names(metadata)[ix], collapse = ", "),
               " files not found for site ", site_nm, "!\n")))
    }

    # read files
    datIn <- lapply(setNames(tbls, nm = sub("\\..*", "", basename(tbls))), read.csv)

    # check data fidelity against template
    err <- 0
    for (i in seq_along(datIn)) {

      # check column names
      err <- checkColNms(metadata, names(datIn)[i], datIn, err)

      # check data types

      # check data values

      # bind to database
      if (err == 0) database[[i]] <- rbind(database[[i]], datIn[[i]])
    }
  }
}
