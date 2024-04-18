# Read ess-dive

read_ess.dive.fx <- function(site_name) {
  p <- file.path(dirname(getwd()), "data", "raw", site_name)
  f <- list.files(p)
  f.ls <- list.files(file.path(p, f[grep("^ess_dive[^.]*$", f)], "data"), full.names = TRUE)
  nms <- basename(f.ls)
  ix <- grep("_dd", nms, invert = TRUE)
  ls <- lapply(f.ls[ix], function(x) 
    read.csv(x, na.strings = c("-9999", "N/A"))[-1, ])
  ls <- lapply(ls, type.convert, as.is = TRUE)
  names(ls) <- sapply(strsplit(nms[ix], split = "\\."), "[[", 1)
  return(ls)
}
