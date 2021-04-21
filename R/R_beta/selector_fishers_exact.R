#' @include helper_functions.R
#-------------------------------------------------------------------------------

#' gof_fisher_exact
#' @title gof_fisher_exact
#' @name gof_fisher_exact
#' @inheritParams gof_transformed
#' @param loss.type One of 'min', 'mean', '0', '1', 'wilcox', and 'ks'.
#' Determines how the loss should be computed from the pvalues of the ks and
#' wilcoxon rank tests.
#' @return loss between 0 and 1. In the best case the loss is 0. 
#' @export gof_fisher_exact
gof_fisher_exact <-
  function(feat,
           yobs,
           tr,
           estimator,
           k = 2,
           verbose = FALSE,
           loss.type = 'min') {
    print("Think about how this could be done when the propensity score is not
          constant.")
    catch_error(feat, yobs, tr, k)
    # ----------------------------------------------------------------------------
    # Compute the CATE estimates using a k fold CV
    cate_est <-
      compute_CATE_estimates(feat, yobs, tr, estimator, k, verbose)
    
    # ----------------------------------------------------------------------------
    # Compute the Y(0) and Y(1)
    
    Y0 <- yobs - ifelse(tr == 1, cate_est, 0)
    Y1 <- yobs + ifelse(tr == 0, cate_est, 0)
    
    pw0 <- wilcox.test(x = Y0[tr == 1], y = Y0[tr == 0])$p.value
    pk0 <- ks.test(x = Y0[tr == 1], y = Y0[tr == 0])$p.value
    pw1 <- wilcox.test(x = Y1[tr == 1], y = Y1[tr == 0])$p.value
    pk1 <- ks.test(x = Y1[tr == 1], y = Y1[tr == 0])$p.value
    
    # ----------------------------------------------------------------------------
    if (loss.type == 'min') {
      loss <- -min(pw0, pk0, pw1, pk1)
    } else if (loss.type == 'mean') {
      loss <- -mean(pw0, pk0, pw1, pk1)
    } else if (loss.type == '0') {
      loss <- -min(pw0, pk0)
    } else if (loss.type == '1') {
      loss <- -min(pw1, pk1)
    } else if (loss.type == 'wilcox') {
      loss <- -min(pw0 , pw1)
    } else if (loss.type == 'ks') {
      loss <- -mean(pk0, pk1)
    }
    return(1 + loss)
  }
