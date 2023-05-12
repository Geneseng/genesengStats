globalVariables(
  c(
    ".",
    "bart",
    "biomarker",
    "cochran",
    "friedman",
    "ftest",
    "glm_roc",
    "group2",
    "kw",
    "mcnemar",
    "mwtest",
    "pairedttest",
    "repeated_anova",
    "Shapiro's test",
    "subject",
    "ttest",
    "value",
    "variable",
    "Var1",
    "Var2",
    "X1",
    "X2",
    "welchtest",
    "wilcoxtest",
    "mann",
    "welch"
    )
)

#' Compute Spearman coefficient and confident intervals at 95\%
#'
#' @param var1 categorical variable
#' @param var2 continuous variable
#'
#' @export
spearmanCI <- function(var1, var2){
  
  coeff <- cor(var1, var2, method = "spearman",  use = "pairwise.complete.obs")
  n <- min(length(var1[!is.na(var1)]), length(var2[!is.na(var2)]))
  stderr <- 1.0 / sqrt(n - 3)
  delta <- 1.96 * stderr
  
  lower <- tanh(atanh(coeff) - delta)
  upper <- tanh(atanh(coeff) + delta)
  
  df <- cbind.data.frame(coeff, lower, upper)
  df1 <- apply(df, 2, round, digits = 3)
  
  return(df1)
  
}
