make_dispersion_fixture <- function() {
  set.seed(42)
  ids <- paste0("s", 1:24)
  counts <- matrix(rpois(20 * 24, 20), 20, dimnames = list(paste0("t", 1:20), ids))
  counts[1:4, 13:24] <- counts[1:4, 13:24] + 40
  map <- data.frame(group = rep(c("A", "B", "C"), each = 8),
                    site = rep(rep(c("X", "Y"), each = 4), 3),
                    block = rep(1:8, 3), covar = rep(c("F", "M"), 12), row.names = ids)
  phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE), phyloseq::sample_data(map))
}

test_that("PERMDISP uses full-space median distances and preserves RNG", {
  ps <- make_dispersion_fixture()
  d <- phyloseq::distance(ps, method = "bray")
  group <- factor(phyloseq::sample_data(ps)$group)
  expected_fit <- vegan::betadisper(d, group, type = "median", bias.adjust = TRUE)
  set.seed(1)
  expected_test <- vegan::permutest(expected_fit, permutations = 999)
  set.seed(73)
  seed <- .Random.seed
  got <- Gotools:::.gotools_permdisp(d, group)
  expect_identical(.Random.seed, seed)
  expect_equal(got$distance, expected_fit$distances)
  expect_equal(got$F, expected_test$tab[1, "F"])
  expect_equal(got$p, expected_test$tab[1, "Pr(>F)"])
})

test_that("Go_bdivPM saves companion PDFs across metrics, facets and combinations", {
  root <- tempfile("permdisp-")
  dir.create(root)
  old_wd <- getwd()
  on.exit({setwd(old_wd); unlink(root, recursive = TRUE)}, add = TRUE)
  setwd(root)
  ps <- make_dispersion_fixture()
  args <- list(psIN = ps, cate.vars = "group", project = "Disp", orders = c("A", "B", "C", "X", "Y"),
               distance_metrics = c("bray", "euclidean"), ellipse = FALSE, height = 4, width = 7,
               mycols = c(A = "#0072B2", B = "#D55E00", C = "#009E73"), addnumber = FALSE)
  suppressWarnings(do.call(Go_bdivPM, args))
  pdf_dir <- sprintf("Disp_%s/pdf", format(Sys.Date(), "%y%m%d"))
  files <- list.files(pdf_dir, full.names = TRUE)
  pm <- files[grepl(".PM.", basename(files), fixed = TRUE)]
  companion <- file.path(pdf_dir, sub(".PM.", ".PERMDISP.", basename(pm), fixed = TRUE))
  expect_length(pm, 1)
  expect_true(file.exists(companion))
  expect_true(all(file.info(files)$size > 1000))
  if (requireNamespace("pdftools", quietly = TRUE)) {
    text <- pdftools::pdf_text(companion)
    expect_match(paste(text, collapse = ""), "PERMDISP")
    expect_match(paste(text, collapse = ""), "spatial median")
  }
  args$name <- "facet"
  args$facet <- "site"
  args$cate.conf <- "covar"
  args$strata_var <- "block"
  args$combination <- 2
  suppressWarnings(do.call(Go_bdivPM, args))
  facet_file <- list.files(pdf_dir, pattern = "PERMDISP.*facet", full.names = TRUE)
  expect_length(facet_file, 1)
  if (requireNamespace("pdftools", quietly = TRUE)) {
    text <- paste(pdftools::pdf_text(facet_file), collapse = "")
    expect_match(text, "covariates not fitted")
    expect_match(text, "PERMDISP")
  }
  args$name <- "patchwork"
  args$facet <- args$cate.conf <- args$strata_var <- args$combination <- NULL
  args$patchwork <- TRUE
  plots <- suppressWarnings(do.call(Go_bdivPM, args))
  expect_length(plots, 2)
  expect_length(list.files(pdf_dir, pattern = "PERMDISP.*patchwork"), 1)
  expect_length(list.files(pdf_dir, pattern = "\\.PM\\..*patchwork"), 0)
  args$name <- "nostats"
  args$statistics <- FALSE
  suppressWarnings(do.call(Go_bdivPM, args))
  expect_length(list.files(pdf_dir, pattern = "PERMDISP.*nostats"), 0)
})

test_that("small groups skip dispersion without breaking ordination", {
  root <- tempfile("permdisp-small-")
  dir.create(root)
  old_wd <- getwd()
  on.exit({setwd(old_wd); unlink(root, recursive = TRUE)}, add = TRUE)
  setwd(root)
  ps <- make_dispersion_fixture()
  ps <- phyloseq::prune_samples(c("s1", paste0("s", 9:16)), ps)
  expect_warning(plots <- withCallingHandlers(
    Go_bdivPM(ps, "group", "Small", c("A", "B"), "euclidean",
              ellipse = FALSE, patchwork = TRUE, addnumber = FALSE),
    warning = function(w) {
      if (grepl("attempt to set", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
    }),
                 "PERMDISP skipped")
  expect_length(plots, 1)
  expect_length(list.files(sprintf("Small_%s/pdf", format(Sys.Date(), "%y%m%d"))), 0)
})
