
# Utility functions

#' Package load hook
#'
#' Initialize package-level options and internal state when the SWEDDIE package
#' is loaded. This function sets default logging behavior and prepares internal
#' environments required by other package functions.
#'
#' This function is called automatically by R when the package namespace is
#' loaded and should not be invoked directly.
#'
#' @param libname character; path to the library directory.
#' @param pkgname character; name of the package.
#'
#' @return NULL (called for side effects only).
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .sweddie_log_opts <- new.env(parent = emptyenv())
  
  .sweddie_log_opts$verbose <- TRUE
  .sweddie_log_opts$append  <- TRUE
  .sweddie_log_opts$file    <- ""
  
  ## assign into package namespace
  assign(".sweddie_log_opts", .sweddie_log_opts,
         envir = parent.env(environment()))
}

#' Configure SWEDDIE logging behavior
#'
#' @param verbose logical; print messages to console?
#' @param file character; path to log file ("" = console only)
#' @param append logical; append to existing log file?
#'
#' @export
sweddie_set_logging <- function(
    verbose = NULL,
    file = NULL,
    append = NULL
) {
  if (!is.null(verbose))
    .sweddie_log_opts$verbose <- isTRUE(verbose)
  
  if (!is.null(file))
    .sweddie_log_opts$file <- as.character(file)
  
  if (!is.null(append))
    .sweddie_log_opts$append <- isTRUE(append)
  
  invisible(.sweddie_log_opts)
}

#' Package-level logging utility
#'
#' Write formatted messages to the console and/or a log file according to the
#' current SWEDDIE logging configuration. Logging behavior (verbosity, output
#' destination, and append mode) is controlled globally via
#' \code{\link{sweddie_set_logging}} and applied consistently across all package
#' functions.
#'
#' This function is intended for internal use only and should not be called
#' directly by users.
#'
#' @param ... Character vectors or objects coercible to character, passed to
#'   \code{\link[base]{cat}}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @seealso \code{\link{sweddie_set_logging}}
#'
#' @examples
#' \dontrun{
#' vcat("Starting data ingestion")
#' vcat("Processed", n, "files")
#' }
vcat <- function(...) {
  opts <- .sweddie_log_opts
  
  if (!isTRUE(opts$verbose))
    return(invisible(NULL))
  
  cat(
    ...,
    "\n",
    file   = opts$file,
    append = opts$append
  )
}

#' Add carriage return to .csv files
#'
#' @param dataDir directory containing .csv files
#' @description Adds a carriage return (end of line character) to .csv files. Note that this is recursive by default. This behavior can be changed by setting "recursive = FALSE".
#' @return modified .csv files
#' @keywords internal
crAdd <- function(dataDir) {
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
#' @note This is typically called only from \code{\link{checkColNms}}.
#' @keywords internal
# checks column names
checkColNms <- function(tableName, datIn, err, file = outfile) {
  
  readCore <- function(coreDir, return = "template") {
    # get dir name as needed
    if (missing(coreDir)) {
      coreDir <- "../data/sweddie/metadata/core"
    }
    # create template list
    ls <- list.files(coreDir, full.names = TRUE)
    ix.dd <- grep("_dd", ls)
    ls.dd <- ls[ix.dd]
    ls.tm <- ls[-ix.dd]
    dd <- lapply(
      setNames(ls.dd, nm = sapply(strsplit(basename(ls.dd), "_"), "[[", 1)),
      read.csv)
    if (return == "template") {
      lapply(dd, function(x) {
        setNames(data.frame(matrix(ncol = nrow(x), nrow = 0)), x[ , 1])
      })
    } else {
      return(dd)
    }
  }
  
  vcat("\t", tableName, "\n", append = TRUE, file = outfile)
  
  # check for template
  metadata <- readCore(return = "dd")
  
  # check that required data files are present (experiment, site, plot)
  if (!any(grepl(tableName, names(datIn)))) {
    err <- err + 1
    vcat("\t", tableName, "is missing in metadata/siteData directory\n", append = TRUE, file = outfile)
  }
  
  # check that required columns are present
  req <- metadata[[tableName]][which(metadata[[tableName]][["req"]] == "yes"), 1]
  for (i in seq_along(req)) {
    if (!any(grepl(req[i], names(datIn[[tableName]])))) {
      err <- err + 1
      vcat("\t\t", "Required column", req[i], "is missing in", tableName, "table\n", append = TRUE, file = outfile)
    }
  }
  
  # check for missing cols in data
  miss <- setdiff(metadata[[tableName]][["col_name"]], colnames(datIn[[tableName]]))
  miss <- miss[!(miss %in% req)]
  
  # check for extra cols in data
  xtra <- setdiff(colnames(datIn[[tableName]]), metadata[[tableName]][["col_name"]])
  
  if (length(miss) > 0 ) {
    vcat("\t\t Non-required columns missing from data:", miss, "\n", append = TRUE, file = outfile)
  }
  if (length(xtra > 0)) {
    xtras <- sapply(seq_along(xtra), function(i) {
      ix <- which(colnames(datIn[[tableName]]) == xtra[i])
      paste0(xtra[i], " (", ix, ")")
    }) 
    vcat("\t\t Extra columns:", xtras, "\n", append = TRUE, file = outfile)
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
read_meta <- function(metaDir, return = "template") {
  # get dir name as needed
  if (missing(metaDir)) {
    metaDir <- "../data/sweddie/metadata/core"
  }
  # create template list
  ls <- list.files(metaDir, full.names = TRUE)
  ls <- ls[grep("_dd", ls)]
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

#' Read template files
#'
#' @param exp_names list of directories containing SWEDDIE data
#' @return A list with core template tables or core meta data
#' @keywords internal
#' @description Wrapper function to run fx \code{\link{read_meta}} on all files in list of directories specified by "exp_names"  
#' @details Each directory in the list 'exp_names' must contain the directories named "input_data" and "meta", per \code{\link{read_meta}} specification
getFLMD_DD <- function(exp_names = NULL) {
  if (missing(exp_names)) {
    exp_names <- list.files(path.expand("~/eco-warm/data/experiments")) 
  }
  flmd_dd.ls <- lapply(
    lapply(exp_names, read_meta, verbose = TRUE), function(x) {
      if (length(x$flmd) == 0) {
        NULL
      } else {
        x
      } 
    })
  names(flmd_dd.ls) <- exp_names
  Filter(Negate(is.null), flmd_dd.ls)
}

#' Get csv names
#'
#' @param path.csv path to csv file
#' @return names of csv file columns
#' @keywords internal
#' @description allows name checking without loading whole file
getCSVnms <- function (path.csv) {
  names(read.csv(path.csv, nrows = 1, check.names = FALSE, strip.white = TRUE))
}

#' Check observation frequency
#' 
#' @param df data frame
#' @param dateName name of column containing observation dates
#' @param pltName name of column with plot identifiers
#' @param depth are data depth resolved? If so, data must contain a column called "depth"
#' @return list
#' @keywords internal
#' @description returns list of median interval in seconds, days, inferred sampling frequency, etc.
ts_freq <- function(df, dateName = "Date", pltName = "plt_name", repName = NULL, depth = FALSE) {
  
  # validate inputs
  if (!dateName %in% names(df)) stop("Column '", dateName, "' not found in data frame.")
  if (!pltName %in% names(df)) stop("Column '", pltName, "' not found in data frame.")
  if (!is.data.frame(df)) stop("Input must be a data frame.")
  
  # filter to single plot
  p1 <- df[[pltName]][1]
  df.f <- df[df[[pltName]] == p1, , drop = FALSE]
  
  # filter to single depth (as needed)
  if (depth && "depth" %in% names(df)) {
    min_depth <- suppressWarnings(min(df.f$depth, na.rm = TRUE))
    df.f <- df.f[df.f$depth == min_depth, , drop = FALSE]
  }
  
  # filter to rep (as needed)
  if (!is.null(repName)) {
    r1 <- df[[repName]][1]
    df.f <- df.f[df.f[[repName]] == r1, , drop = FALSE]
  }
  
  # Parse and clean date column
  # Convert to POSIXct or Date robustly
  dates_raw <- df.f[[dateName]]
  
  # Handle different possible formats automatically
  if (!inherits(dates_raw, c("Date", "POSIXct", "POSIXt"))) {
    dates <- suppressWarnings(lubridate::parse_date_time(
      dates_raw,
      orders = c("Ymd HMS", "Ymd HM", "Ymd H", "Ymd", "mdY", "dmy", "ymd"),
      tz = "UTC"
    ))
  } else {
    dates <- as.POSIXct(dates_raw)
  }
  
  # Drop invalid or NA dates
  dates <- sort(dates[!is.na(dates)])
  
  if (length(dates) < 2) {
    warning("Not enough valid dates to estimate frequency.")
  }
  
  # Compute intervals
  diffs <- diff(dates)
  dt <- as.numeric(stats::median(diffs, na.rm = TRUE), units = "secs")
  
  # --- 5. Map to human-readable frequency -----------------------------------
  freq <- dplyr::case_when(
    dt < 60 ~ "sub-minute",
    dt < 3600 ~ "minutely",
    dt < 86400 ~ "hourly",
    dt < 86400 * 7 ~ "daily",
    dt < 86400 * 31 ~ "weekly",
    dt < 86400 * 365 ~ "monthly",
    TRUE ~ "yearly"
  )
  
  # Return clean summary
  list(
    median_interval_days = round(dt / 86400, 3),
    inferred_frequency = freq,
    date_range = paste0(format(dates[1]), " to ", format(dates[length(dates)])),
    n_obs = length(dates)
  )
}

#' Build core SWEDDIE database
#' 
#' @param DIR parent directory in which SWEDDIE 'database' directory is stored
#' @param write_report should report of build be written to a file?
#' @param verbose should output be printed to console?
#' @param append call for internal utility function 'vcat'
#' @return list
#' @keywords internal
#' @description returns SWEDDIE core database object with meta, site, and plot tables
coreData.fx <- function(DIR = "../data/sweddie", write_report = TRUE, verbose = TRUE) {
  
  # Constants
  DB_DIR <- "database"
  S_DIR <- "siteData"
  LIST_FILE <- "coreData.rda"
  TIMESTAMP <- format(Sys.time(), "%y%m%d-%H%M") 
  
  # configure logging for this run
  .sweddie_log_opts$verbose <- verbose
  .sweddie_log_opts$append  <- TRUE
  
  # Set output file
  outfile <- ""
  if (write_report) {
    outfile <- file.path(DIR, DB_DIR, paste0("logs/coreLog", "_", TIMESTAMP, ".txt"))
    invisible(file.create(outfile))
    .sweddie_log_opts$file <- outfile
  } else {
    .sweddie_log_opts$file <- ""  # console only
  }
  
  # Start writing in the output file
  vcat("SWEDDIE Compilation Log \n",
       "\n", as.character(Sys.time()),
       "\n", rep("-", 15), "\n", file = outfile, verbose = verbose)
  
  vcat("\n\nCompiling data files in", S_DIR, "\n", rep("-", 30), "\n", file = outfile, append = append, verbose = verbose)
  
  data_dirs <- list.dirs(file.path(DIR, S_DIR), full.names = TRUE, recursive = FALSE)
  if (!length(data_dirs)) {
    vcat("No data directories found!\n", file = outfile, append = append, verbose = verbose)
    return(NULL)
  }
  
  # ensure EOL carriage return present
  invisible(crAdd(file.path(DIR, S_DIR)))
  
  vcat("Compiling and checking core data...\n\n", file = outfile, append = append, verbose = verbose)
  #   pb <- txtProgressBar(min = 0, max = length(data_dirs), style = 3)
  # }
  # 
  # check if previous database object exists in database directory, and only update file if new data exisit
  if (file.exists(file.path(DIR, DB_DIR, LIST_FILE))) {
    
    # load existing database
    load(file.path(dataset_directory, DB_DIR, LIST_FILE)) # obj "coreDat"
    
    # convert to character and coerce to list of data frames
    coreDat_chr <- lapplydf(lapply(coreData, function(x) lapply(x, as.character)))
    
    # remove old version
    rm(coreDat)
    
    # Split each table by entry_name
    coreDat_old <- lapply(coreDat_chr, function(x) split(x, x$exp_name))
  } else {
    database <- setNames(vector(mode = "list", length = length(data_dirs)), nm = basename(data_dirs))
  }
  
  # compile new templates and check against existing data
  for (d in seq_along(database)) {
    cat(names(database[d]), "\n")
    # get expName
    expName <- names(database[d])
    
    vcat("\n", expName, "\n", file = outfile, append = append, verbose = verbose)
    
    # get tables
    tbls <- list.files(data_dirs[d], full.names = TRUE)
    
    # read files
    datIn <- lapply(
      setNames(tbls, nm = sub("\\..*", "", basename(tbls))), 
      function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE, strip.white = TRUE)
        
        # Trim leading/trailing spaces in all character columns
        df[] <- lapply(df, function(x) {
          if (is.character(x)) trimws(x) else x
        })
        
        # Remove rows where all values are NA or empty string
        df <- df[rowSums(is.na(df) | df == "") != ncol(df), ]
        
        df
      }
    )
    
    # check data fidelity against template
    err <- 0
    for (i in seq_along(datIn)) {
      
      # check column names
      err <- checkColNms(names(datIn)[i], datIn, err, file = outfile)
      
      # check data types
      
      # check data values
    }
    
    # bind to database
    if (err == 0) database[[d]] <- datIn
  }
  database
}
