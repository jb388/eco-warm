# unit tests

# 1) Logging respects verbose = FALSE
test_that("vcat is silent when verbose = FALSE", {
  sweddie_set_logging(verbose = FALSE)
  
  out <- capture.output(vcat("Hello world"))
  
  expect_length(out, 0)
})

# 2) Logging writes to console when verbose = TRUE
test_that("vcat prints when verbose = TRUE", {
  sweddie_set_logging(verbose = TRUE, file = "")
  
  out <- capture.output(vcat("Hello"))
  
  expect_true(any(grepl("Hello", out)))
})

# 3) Logging writes to file
test_that("vcat writes to file", {
  withr::with_tempfile("log", {
    sweddie_set_logging(
      verbose = TRUE,
      file = log,
      append = TRUE
    )
    
    vcat("File test")
    
    expect_true(file.exists(log))
    expect_true(any(grepl("File test", readLines(log))))
  })
})

# 4) append works
test_that("vcat appends to existing file", {
  withr::with_tempfile("log", {
    sweddie_set_logging(verbose = TRUE, file = log, append = TRUE)
    
    vcat("First")
    vcat("Second")
    
    lines <- readLines(log)
    expect_true(any(grepl("First", lines)))
    expect_true(any(grepl("Second", lines)))
  })
})

# coreData.fx logging
test_that("coreData.fx sets logging context", {
  withr::with_tempdir({
    expect_silent(
      coreData.fx(DIR = tempdir(), write_report = FALSE, verbose = FALSE)
    )
  })
})