#' gclr.test: Test arbitrary weighted taxa groups (GCLR inference)
#'
#' @description
#' This function performs statistical hypothesis testing on taxa groups
#' (with arbitrary weighted components) for GCLR inference. It identifies taxa
#' groups that are significant in relative abundance (relative to the user-defined
#' reference taxa group).
#'
#' @usage
#' gclr.test(model, variables, w1, w2)
#'
#' @param model The object fitted using function \code{McDisCov}.
#' @param variables Vector containing variables of interest.
#'   For example: \code{variables = c("binDisease","cont")}.
#'   For categorical variables, make sure the names align with the column names of \code{model$meta_data}.
#' @param w1 Vector; weights of the taxa in the target taxa group. It will be scaled to unit sum automatically.
#' @param w2 Vector; weights of the taxa in the reference taxa group. It will be scaled to unit sum automatically.
#'
#' @return
#' A data frame with columns:
#' \itemize{
#' \item \code{covariate_name}: Variables of interest.
#' \item \code{estimate}: Effect size estimate in log scale for the target group relative to the reference group.
#' \item \code{SE}: Standard error of the estimate.
#' \item \code{TestStat}: Test statistic (chi-squared, 1 df).
#' \item \code{p_value}: p-value.
#' }
#'
#' @examples
#' \dontrun{
#' # CLR
#' mcdiscov_clr <- function(model, variables) {
#'   res <- data.frame()
#'   for (i in seq(ncol(model$taxa_data))) {
#'     w1 <- rep(0, ncol(model$taxa_data))
#'     w2 <- rep(1, ncol(model$taxa_data))
#'     w1[i] <- 1
#'     res_i <- gclr.test(model, variables, w1 = w1, w2 = w2)
#'     res <- rbind(res, res_i)
#'   }
#'   rownames(res) <- colnames(model$taxa_data)
#'   return(res)
#' }
#'
#' # ALR
#' mcdiscov_alr <- function(model, variables, reference = 1) {
#'   res <- data.frame()
#'   for (i in seq(ncol(model$taxa_data))) {
#'     w1 <- rep(0, ncol(model$taxa_data))
#'     w2 <- rep(0, ncol(model$taxa_data))
#'     w1[i] <- 1
#'     w2[reference] <- 1
#'     res_i <- gclr.test(model, variables, w1 = w1, w2 = w2)
#'     res <- rbind(res, res_i)
#'   }
#'   rownames(res) <- colnames(model$taxa_data)
#'   return(res)
#' }
#' }
#'
#' @author Jiahao Wang (\email{10225000478@stu.ecnu.edu.cn})
#' @export
gclr.test <- function(model, variables = NULL, w1, w2) {

  # prep
  analysis_result <- model$numerical_res
  X <- model$meta_data
  p <- length(w1)
  stopifnot(length(w2) == p)
  q <- ncol(X)

  alpha_taxmaj <- analysis_result$estimated_theta
  Sigma_all    <- analysis_result$theta_cov_matrix
  stopifnot(length(alpha_taxmaj) == p*q, all(dim(Sigma_all) == c(p*q, p*q)))

  # start
  if (is.null(variables)) {
    variables_idx <- 1:q
  } else if (is.character(variables)) {
    stopifnot(!is.null(colnames(X)))
    variables_idx <- match(variables, colnames(X))
    if (any(is.na(variables_idx))) {
      miss <- variables[is.na(variables_idx)]
      stop(sprintf("Variables not found: %s", paste(miss, collapse = ", ")))
    }
  } else {
    variables_idx <- as.integer(variables)
    stopifnot(all(variables_idx >= 1 & variables_idx <= q))
  }

  # weighting
  norm_group <- function(w) {
    s <- sum(w)
    if (s > 0) w / s else w
  }
  w1n <- norm_group(w1)
  w2n <- norm_group(w2)
  c_vec <- w1n - w2n

  # calc
  out_list <- lapply(variables_idx, function(k) {
    idx_k   <- seq(k, p*q, by = q)
    alpha_k <- alpha_taxmaj[idx_k]
    Sigma_k <- Sigma_all[idx_k, idx_k] # p×p

    est <- sum(c_vec * alpha_k)
    se2 <- as.numeric(t(c_vec) %*% Sigma_k %*% c_vec)
    se  <- sqrt(max(se2, 0))
    chi2 <- if (se > 0) (est / se)^2 else 0
    pval <- pchisq(chi2, df = 1, lower.tail = FALSE)

    data.frame(
      covariate_name  = if (!is.null(colnames(X))) colnames(X)[k] else NA_character_,
      estimate = est,
      SE       = se,
      TestStat = chi2,
      p_value  = pval,
      row.names = NULL
    )
  })

  do.call(rbind, out_list)
}
