# Use the full distance space, not distances measured on two plotted axes.
.gotools_permdisp <- function(d, group, permutations = 999) {
  seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (is.null(seed)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", seed, envir = .GlobalEnv))
  set.seed(1)
  fit <- vegan::betadisper(d, group, type = "median", bias.adjust = TRUE)
  test <- vegan::permutest(fit, permutations = permutations)
  list(distance = fit$distances, F = test$tab[1, "F"], p = test$tab[1, "Pr(>F)"])
}
