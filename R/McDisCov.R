#' McDisCov: Outcome-linked Microbiome Discovery via Compositional lOgit-normal modeling under Varying support
#'
#' @description
#' The function identifies taxa whose absolute abundances (actual counts in a unit volume of an ecosystem)
#' and relative abundances (actual counts relative to the auto-selected reference taxa group) are associated
#' with the variable(s) of interest. It supports multiple continuous and categorical covariates for association
#' identification. For additional identification of arbitrarily weighted taxa groups (such as ALR, CLR, ILR,
#' and balance), refer to the function \code{gclr.test}.
#'
#' @usage
#' library(McDisCov)
#' McDisCov(
#'   taxa_data,
#'   meta_data,
#'   variables,
#'   imputation = FALSE,
#'   zero_impute = 0.5,
#'   zero_ratio_threshold = 1,
#'   alpha = 0.05,
#'   Wc = 0.7,
#'   nonzero_cut = 3,
#'   lib_cut = 0,
#'   cov_shrinkage = "diag",
#'   lambda_seq = seq(0, 1, by = 0.1)
#' )
#'
#' @param taxa_data A data frame or matrix representing the observed count table.
#'   Rows correspond to samples and columns correspond to taxa.
#' @param meta_data A data frame of covariates. The rows of \code{meta_data} correspond
#'   to the rows of \code{taxa_data}.
#' @param variables A character vector containing covariates of interest and potential
#'   confounders. For example, \code{variables = c("bin", "cont")}. Note that grouping
#'   variables should be converted into factors.
#' @param imputation Logical; \code{TRUE} or \code{FALSE}. Default is \code{FALSE}.
#' @param selected_taxa Integer indices of reference taxa to always keep or impute;
#'   default is \code{NULL} (automatic selection).
#' @param zero_impute A positive real value. Default is \code{0.5}. If
#'   \code{imputation = TRUE}, zeros are replaced with this pseudo-count.
#' @param zero_ratio_threshold A real number in the interval \eqn{[0,1]}. Default is \code{1}.
#'   When \code{imputation = TRUE}, taxa are sorted by zero proportion, and those with zero ratios
#'   less than or equal to \code{zero_ratio_threshold} will be imputed with \code{zero_impute}.
#' @param alpha A real value between 0 and 1; the significance level for differential
#'   relative abundance tests. Default is \code{0.05}.
#' @param Wc A real value between 0 and 1; the cutoff used to identify taxa with weak and
#'   strong associations. Default is \code{0.7}.
#' @param nonzero_cut A non-negative integer. Taxa with more than \code{nonzero_cut}
#'   nonzero observations per group are kept (a continuous variable is treated as one group).
#'   Default is \code{3}.
#' @param lib_cut A non-negative real value. After taxa filtering, samples with total read counts
#'   greater than \code{lib_cut} are retained. Default is \code{0}.
#' @param cov_shrinkage Character string specifying the shrinkage target for covariance estimation.
#'   Default is \code{"diag"}. Options include \code{"diag"}, \code{"mean"}, and \code{"identity"}.
#'   See the paper for more details.
#' @param lambda_seq A numeric vector specifying the shrinkage strengths to be tested.
#'   Default is \code{seq(0, 1, by = 0.1)}.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{taxa_data}: The taxa table used for analysis.
#'   \item \code{meta_data}: The metadata used for analysis.
#'   \item \code{AA}: A named list (one element per tested variable). Each element is a data frame with columns:
#'     \describe{
#'       \item{Taxa}{Taxon name.}
#'       \item{effect_size}{Estimated log-scale effect size with respect to absolute abundance.}
#'       \item{W}{Total number of significant pairwise differences for that taxon.}
#'       \item{DA_group}{Association strength per taxon: \code{0 = not}, \code{1 = weak}, \code{2 = strong}.}
#'     }
#'   \item \code{RA}: A named list (one element per tested variable). Each element is a data frame with columns:
#'     \describe{
#'       \item{Taxa}{Taxon name.}
#'       \item{Estimate}{Per-taxon log-scale effect size estimate relative to the selected reference taxa.}
#'       \item{StdError}{Standard error of the estimate.}
#'       \item{TestStatistic}{Chi-squared test statistic with 1 degree of freedom.}
#'       \item{p_val}{P-value.}
#'       \item{DA}{Differential-abundance flag relative to the reference taxa: \code{0 = not associated}, \code{1 = associated}.}
#'     }
#'   \item \code{selected_taxa}: The automatically selected reference taxa.
#'   \item \code{numerical_res}: A list containing numerical outputs:
#'     \code{estimated_theta}, \code{estimated_sigma_c}, \code{theta_cov_matrix}, \code{pw_est}, and \code{pw_sig}.
#' }
#'
#' @author
#' Jiahao Wang (\email{10225000478@stu.ecnu.edu.cn})
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' Z <- matrix(rpois(100 * 20, 5), 100, 20)
#' meta <- data.frame(group = factor(sample(c("A", "B"), 100, TRUE)))
#' res <- McDisCov(Z, meta, variables = "group")
#' str(res$AA$group)
#' }
#'
#' @export


#-----------------------------------------------------------------------------------------------------------------------------------------------
McDisCov <- function(
    taxa_data,                      # samples x taxa counts (matrix or data.frame)
    meta_data,                      # samples x covariates (data.frame; same rows as taxa_data)
    variables,                      # character vector of covariate names to include & test
    imputation = FALSE,
    selected_taxa = NULL,
    zero_ratio_threshold = 1,
    zero_impute = 0.5,
    Wc = 0.7,
    alpha = 0.05,
    nonzero_cut = 3,
    lib_cut = 0,
    cov_shrinkage = "diag",
    lambda_seq = seq(0, 1, by = 0.1)
){
  ## --------------------------- 1) Data checks & design matrix ---------------------------
  message("Data preprocessing...")

  keep_list <- data_preprocessing(taxa_data = taxa_data, meta_data = meta_data,
                                  nonzero_cut = nonzero_cut, lib_cut = lib_cut)

  Z <- taxa_data[keep_list$sample.keep, keep_list$taxa.keep, drop = FALSE]
  X <- meta_data[keep_list$sample.keep, , drop = FALSE]

  # coerce Z (counts) to numeric matrix
  Z <- as.matrix(Z)
  mode(Z) <- "numeric"

  # row alignment
  if (nrow(Z) != nrow(meta_data)) {
    stop("Rows of taxa_data and meta_data must match (same samples, same order).")
  }

  # check variables
  if (!all(variables %in% colnames(meta_data))) {
    miss <- setdiff(variables, colnames(meta_data))
    stop("variables not found in meta_data: ", paste(miss, collapse = ", "))
  }

  message("Building design matrix X from meta_data[variables] ...")
  X_raw <- meta_data[, variables, drop = FALSE]

  # Expand factors to 0/1 dummies (no intercept)
  if (any(!vapply(X_raw, is.numeric, logical(1)))) {
    X_vars <- model.matrix(~ . , data = X_raw)[, -1, drop = FALSE]  # dummy expansion, no intercept
  } else {
    X_vars <- as.matrix(X_raw)
  }
  storage.mode(X_vars) <- "numeric"

  # Add intercept as the LAST column
  X <- cbind(X_vars, Intercept = 1)

  ## --------------------------- 2) Zero handling / reference selection -------------------
  message("Selecting reference / zero handling ...")
  zero_ratio <- colSums(Z == 0) / nrow(Z)

  if (!is.null(selected_taxa)) {
    for (j in selected_taxa) Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
  } else {
    if (imputation) {
      selected_taxa <- which(zero_ratio <= zero_ratio_threshold)
      if (length(selected_taxa) == 0) {
        stop("No taxon selected. Choose another zero_ratio_threshold.")
      }
      for (j in selected_taxa) Z[, j] <- ifelse(Z[, j] == 0, zero_impute, Z[, j])
    } else {
      if (length(which(zero_ratio == 0)) == 0) {
        # pick a taxon with fewest zeros, drop samples where it is zero,
        # and ensure no covariate column becomes all-0 or all-1
        sorted_indices <- order(zero_ratio)
        i <- 1
        repeat {
          least0_index <- sorted_indices[i]
          remove_idx <- which(Z[, least0_index] == 0)
          Z_new <- Z[-remove_idx, , drop = FALSE]
          X_new <- X[-remove_idx, , drop = FALSE]

          X_tmp <- X_new[, -ncol(X_new), drop = FALSE]  # drop intercept
          is_all_0_or_1 <- apply(X_tmp, 2, function(col) all(col == 0) || all(col == 1))

          if (any(is_all_0_or_1)) {
            i <- i + 1
            if (i > length(sorted_indices)) {
              stop("No suitable taxon; all possibilities make some X column all 0 or all 1.")
            }
          } else {
            break
          }
        }
        selected_taxa <- sorted_indices[i]
        Z <- Z_new
        X <- X_new
      } else {
        selected_taxa <- which(zero_ratio == 0)
      }
    }
  }

  ## --------------------------- 3) Transform & core analysis -----------------------------
  message("GCLR transforming ...")
  Z_trans <- GCLR_transform(Z, selected_taxa, zero_impute)

  analysis_result <- run_analysis_core(Z_trans, X, selected_taxa, cov_shrinkage, lambda_seq)

  ## --------------------------- 4) Wald tests & grouping --------------------------------
  p <- ncol(Z)   # #taxa
  q <- ncol(X)   # #covariates incl. intercept

  which_X_index <- setdiff(seq_len(q), q)   # by default: test all except intercept

  AA <- list()   # <-- new
  RA <- list()   # <-- new

  taxa_names <- colnames(Z)

  for (k in which_X_index) {
    idx_range <- seq(k, p * q, by = q)
    alpha_extraction_index <- (k*p - p + 1):(k*p)

    est_alpha_this <- analysis_result$estimated_alpha[alpha_extraction_index]
    cov_this       <- analysis_result$cov_matrix[idx_range, idx_range]
    std_error      <- sqrt(diag(cov_this))

    ## --- Wald stats ---
    wald_stat <- (est_alpha_this / std_error)^2
    wald_stat <- ifelse(wald_stat == Inf, 0, wald_stat)
    pval      <- 1 - pchisq(wald_stat, df = 1)
    p_star    <- ifelse(pval < alpha, 1, 0)     # 0.05 default

    ## --- Pairwise differences ---
    var_vec  <- diag(cov_this)
    D        <- outer(est_alpha_this, est_alpha_this, "-")
    Var_mat  <- outer(var_vec, var_vec, "+") - 2 * cov_this
    Wmat     <- D / sqrt(Var_mat)
    Pmat     <- 1 - pchisq(Wmat^2, df = 1)

    wald_test_matrix <- (Pmat < alpha) * 1     # 0.05 default
    wald_test_matrix[is.na(wald_test_matrix)] <- 0

    row_sums <- rowSums(wald_test_matrix)   # total significant pairwise differences per taxon

    ## --- (A/B/C) gives group_star = 0/1/2 ---
    hca_count       <- hclust(dist(row_sums), "complete")
    hca_count_index <- cutree(hca_count, k = 2)

    min_group_index <- hca_count_index[which.min(row_sums)]
    max_group_index <- hca_count_index[which.max(row_sums)]
    Group_A  <- which(hca_count_index == min_group_index)
    Group_BC <- which(hca_count_index == max_group_index)
    Group_B  <- Group_BC[which(row_sums[Group_BC] <  (p * Wc))]
    Group_C  <- Group_BC[which(row_sums[Group_BC] >= (p * Wc))]

    group_star <- rep(0, p)  # 0 = not, 1 = weak, 2 = strong
    group_star[Group_B] <- 1
    group_star[Group_C] <- 2

    ## --- Effect size relative to Group A mean (as in your code) ---
    effect_size <- est_alpha_this - mean(est_alpha_this[Group_A])

    var_name <- colnames(X)[k]

    ## ===================== PACK NEW OUTPUTS =====================

    # 1) AA: per-variable summary for association strength
    #    columns: Taxa, effect_size, DA_group (0/1/2), W (row_sums)
    AA[[var_name]] <- data.frame(
      Taxa        = taxa_names,
      effect_size = effect_size,
      DA_group    = as.integer(group_star),
      W           = as.numeric(row_sums),
      row.names   = NULL
    )

    # 2) RA: per-variable regression table
    #    columns: Taxa, Estimate, StdError, TestStatistic, p_val, DA (0/1)
    RA[[var_name]] <- data.frame(
      Taxa          = taxa_names,
      Estimate      = est_alpha_this,
      StdError      = std_error,
      TestStatistic = wald_stat,          # (χ^2 with 1 df). If you prefer Z, use est/std instead.
      p_val         = pval,
      DA            = as.integer(p_star), # 0/1
      row.names     = NULL
    )
  }

  ## --------------------------- 5) Numerical results bundle ------------------------------
  estimated_theta    <- analysis_result$estimated_alpha
  estimated_theta    <- as.vector(t(matrix(estimated_theta, nrow = p, ncol = q))) # rearrange
  estimated_sigma_c  <- analysis_result$estimated_sigma_c
  theta_cov_matrix   <- analysis_result$cov_matrix

  res <- list(
    meta_data      = X,
    taxa_data      = Z,
    AA             = AA,
    RA             = RA,
    selected_taxa  = selected_taxa,
    numerical_res  = list(
      estimated_theta   = estimated_theta,
      estimated_sigma_c = estimated_sigma_c,
      theta_cov_matrix  = theta_cov_matrix
    )
  )

  return(res)
}


run_analysis_core <- function(Z, X, selected_taxa, cov_shrinkage, lambda_seq) {
  X <- as.matrix(X)
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Z)

  nonzero_index <- lapply(1:n, function(i) which(!is.na(Z[i, ])))

  PhiiXi_list <- lapply(1:n, function(i) {
    row_vector <- X[i, , drop = FALSE]
    identity_p <- diag(p)
    kronecker(identity_p, row_vector)[nonzero_index[[i]], , drop = FALSE]
  })

  # initial estimation
  get_residual_matrix <- function(Z, X) {
    n <- nrow(Z)
    p <- ncol(Z)
    resid_mat <- matrix(NA, nrow = n, ncol = p)

    for (j in 1:p) {
      y_j <- Z[, j]
      non_missing <- !is.na(y_j)

      # zero separation
      fit <- lm(y_j[non_missing] ~ X[non_missing, , drop = FALSE])
      resid_mat[non_missing, j] <- residuals(fit)
    }

    return(resid_mat)
  }

  # Loss function (log-likelihood)
  log_lik <- function(Sigma_c) {
    objective_value <- 0

    for (i in 1:n) {
      nonzero_i <- nonzero_index[[i]]
      PhiiZi <- Z[i, nonzero_i]
      PhiiSigma_cPhiiT <- Sigma_c[nonzero_i, nonzero_i]

      Ri <- chol(PhiiSigma_cPhiiT)

      log_det_i <- 2 * sum(log(diag(Ri)))
      v <- forwardsolve(t(Ri), PhiiZi)
      quad <- sum(v^2)

      # Update total objective
      objective_value <- objective_value + (-0.5) * log_det_i + (-0.5) * quad
    }
    objective_value
  }

  # Sigma
  generate_covariance_matrix <- function() {
    delta_mat <- ifelse(is.na(Z), 0, 1)
    residual <- get_residual_matrix(Z, X)

    # Residual Covariance
    Sigma_emp <- matrix(0, p, p)
    count_mat <- matrix(0, p, p)

    for (i in 1:n) {
      delta_i <- delta_mat[i, ]
      r_i <- residual[i, ]
      if (sum(delta_i) == 0) next

      r_i[delta_i == 0] <- 0
      Sigma_emp <- Sigma_emp + tcrossprod(r_i)
      count_mat <- count_mat + tcrossprod(delta_i)
    }

    Sigma_emp <- ifelse(count_mat > 0, Sigma_emp / count_mat, 0)

    # function is_pos
    is_posdef <- function(Sigma, tol = 1e-8) {
      if (!isSymmetric(Sigma)) return(FALSE)
      eigvals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      all(eigvals > tol)
    }

    best_loglik <- -Inf
    best_Sigma <- NULL
    loglik_values <- length(lambda_seq)

    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]

      # shrink to diag
      if (cov_shrinkage == "diag") {
        Tm <- diag(diag(Sigma_emp))
      } else if (cov_shrinkage == "mean") {
        Tm <- diag(mean(diag(Sigma_emp)), p)
      } else if (cov_shrinkage == "identity") {
        Tm <- diag(1, p)
      } else {
        stop("Unknown cov_shrinkage")
      }

      Sigma_c <- (1 - lambda) * Sigma_emp + lambda * Tm

      # Numerical Correction
      Sigma_c <- if (is_posdef(Sigma_c)) Sigma_c else as.matrix(Matrix::nearPD(Sigma_c)$mat)

      ll <- log_lik(Sigma_c)
      loglik_values[i] <- ll

      if (ll > best_loglik) {
        best_loglik <- ll
        best_Sigma <- Sigma_c
      }
    }
    return(list(
      best_lambda = lambda_seq[which.max(loglik_values)],
      best_loglik = best_loglik,
      best_cov = best_Sigma,
      all_loglik = loglik_values,
      lambda_seq = lambda_seq
    ))
  }

  # constrained theta
  s <- rep(0, p)
  s[selected_taxa] <- 1
  q <- ncol(X)
  C <- kronecker(t(s), diag(q))
  solve_constrained_theta <- function(H, b, C) {
    aug_matrix <- rbind(
      cbind(H, t(C)),
      cbind(C, matrix(0, nrow = q, ncol = q))
    )
    rhs <- rbind(b, matrix(0, nrow = q, ncol = 1))
    sol <- qr.solve(aug_matrix, rhs, tol = 1e-20)
    theta_hat <- sol[1:nrow(H), , drop = FALSE]
    return(theta_hat)
  }

  optimize_alternating <- function(max_iter, tol) {
    message("Estimating covariance...")
    Sigma_c <- generate_covariance_matrix()$best_cov
    Sigma_list <- lapply(1:n, function(i) Sigma_c[nonzero_index[[i]], nonzero_index[[i]]])

    message("Estimating effect size...")
    PhiiZi_list <- lapply(1:n, function(i) Z[i, nonzero_index[[i]]])
    result <- fast_half1_half2(PhiiXi_list, PhiiZi_list, Sigma_list)
    half1 <- result$half1
    half2 <- matrix(result$half2, ncol=1)
    theta_hat <- solve_constrained_theta(half1, half2, C)
    alpha_vector <- as.vector(t(matrix(theta_hat, nrow = q, ncol = p)))

    list(alpha_vector = alpha_vector, Sigma_c = Sigma_c)
  }

  result <- optimize_alternating(max_iter, tol)

  message("Conducting hypothesis testing...")
  compute_full_hessian <- function() {
    alpha_hat <- result$alpha_vector
    Sigma_c_hat <- result$Sigma_c

    Sigma_list <- lapply(1:n, function(i) Sigma_c_hat[nonzero_index[[i]], nonzero_index[[i]]])

    # Rcpp
    H_theta_theta <- fast_H_theta_theta(PhiiXi_list, Sigma_list)
    H_total <- -H_theta_theta

    return(H_total)
  }
  H_full <- compute_full_hessian()
  # Numerical Correction
  if (kappa(H_full) > 1e12) {
    H_full <- H_full + diag(1e-8, nrow(H_full))
  }
  JC <- t(C)
  inv_Hessian <- -qr.solve(H_full, tol = 1e-20)

  cov_matrix <- inv_Hessian - inv_Hessian %*% JC %*%
    qr.solve(t(JC) %*% inv_Hessian %*% JC, tol = 1e-20) %*% t(JC) %*% inv_Hessian
  # Numerical Correction
  diag(cov_matrix)[diag(cov_matrix) < 0] <- 0

  estimated_alpha <- result$alpha_vector
  estimated_sigma_c <- result$Sigma_c

  return(list(
    estimated_alpha = estimated_alpha,
    estimated_sigma_c = estimated_sigma_c,
    cov_matrix = cov_matrix
  ))
}

GCLR_transform <- function(matrix1, selected_taxa, zero_impute) {
  U <- ifelse(matrix1 == 0, 0, 1)
  U[, selected_taxa] <- 1
  NA_U <- ifelse(U == 0, NA, 1)
  p <- ncol(matrix1)
  D <- length(selected_taxa)
  matrix2 <- ifelse(matrix1 == 0, zero_impute, matrix1)
  L <- log(matrix2)
  weights <- t(rep(0, p))
  omega <- 100^2 / (1 + 100^2 * D) # 1/D
  weights[selected_taxa] <- omega
  Gamma <- diag(1, p) - rep(1, p) %*% weights
  matrix2 <- L %*% t(Gamma)
  matrix2 <- NA_U * matrix2
  return(matrix2)
}

data_preprocessing <- function(taxa_data = NULL, meta_data = NULL, nonzero_cut = 3, lib_cut = 0) {
  # taxa filtering
  p <- ncol(taxa_data)
  taxa.keep <- rep(TRUE, p)

  for (var in names(meta_data)) {

    if (is.factor(meta_data[[var]])) {
      # -----------------------------
      # categorical (factor)
      # -----------------------------
      levels_var <- levels(meta_data[[var]])
      for (j in seq_len(p)) {
        taxon_col <- taxa_data[, j]
        valid_all_levels <- all(sapply(levels_var, function(lv) {
          idx <- which(meta_data[[var]] == lv)
          sum(taxon_col[idx] > 0) > nonzero_cut
        }))
        if (!valid_all_levels) {
          taxa.keep[j] <- FALSE
        }
      }

    }
    else {
      # -----------------------------
      # continuous
      # -----------------------------
      for (j in seq_len(p)) {
        taxon_col <- taxa_data[, j]
        if (sum(taxon_col > 0) <= nonzero_cut) {
          taxa.keep[j] <- FALSE
        }
      }
    }
  }

  # sample filtering
  sample.keep <- rowSums(taxa_data[ ,taxa.keep]) > lib_cut

  # keep list
  keep_list <- list(taxa.keep = taxa.keep,
                    sample.keep = sample.keep)
  return(keep_list)
}

