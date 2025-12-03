custom_p_format <- function(p) {
  rstatix::p_format(p, accuracy = 0.001, digits = 2, leading.zero = TRUE)
}

run_pgls_features <- function(compdat, features, response = "Compartment", level = 0.95) {
  
  results <- map_df(features, function(feat) {
    # build formula dynamically
    f <- as.formula(paste(feat, "~", response))
    
    # fit model
    model <- caper::pgls(f, data = compdat)
    
    co  <- summary(model)$coefficients
    est <- co[, "Estimate"]
    se  <- co[, 2]
    df  <- nobs(model) - length(est)            # residual df
    crit <- qt(1 - (1 - level)/2, df = df)
    
    tibble(
      feature    = feat,
      term       = rownames(co)[2],
      estimate   = est[2],
      std_error  = se[2],
      statistic  = co[2,3],
      p_value    = co[2,4],
      r_squared  = summary(model)$r.squared,
      lambda     = model$param["lambda"],
      df = df,
      conf_low  = estimate - crit * std_error,
      conf_high = estimate + crit * std_error
    )
  })
  
  return(results)
}


# Phylogenetic LM ---------------------------------------------------------

run_phylolm_features <- function(df, tree, features, response = "Compartment", level = 0.95) {
  
  results <- purrr::map_df(features, function(feat) {
    # build formula dynamically
    f <- as.formula(paste(feat, "~", response))
    
    # fit model
    model <- phylolm::phylolm(f, data = df, phy = tree, model = "lambda")
    
    co  <- summary(model)$coefficients
    est <- co[, "Estimate"]
    se  <- co[, 2]
    df  <- nobs(model) - length(est)
    crit <- qt(1 - (1 - level)/2, df = df)
    
    tibble(
      feature    = feat,
      term       = rownames(co)[2],
      estimate   = est[2],
      std_error  = se[2],
      statistic  = co[2,3],
      p_value    = co[2,4],
      r_squared  = summary(model)$adj.r.squared,
      lambda     = summary(model)$optpar,
      df = df,
      conf_low  = estimate - crit * std_error,
      conf_high = estimate + crit * std_error
    )
  })
  
  return(results)
}

sig_one <- function(tree, dat, feat_col) {
  x <- dat[[feat_col]]
  names(x) <- dat$strain
  
  # drop NAs and tips accordingly
  ok <- !is.na(x)
  if (sum(ok) < 5) return(tibble(feature = feat_col, n = sum(ok),
                                 lambda = NA_real_, p_value = NA_real_))
  tr <- ape::drop.tip(tree, setdiff(tree$tip.label, names(x)[ok]))
  xv <- x[ok]
  
  # try both metrics
  lam <- try(phytools::phylosig(tr, xv, method = "lambda", test = TRUE), silent = TRUE)

  tibble(
    feature   = feat_col,
    n         = length(xv),
    lambda    = if (inherits(lam, "try-error")) NA_real_ else lam$lambda,
    p_value = if (inherits(lam, "try-error")) NA_real_ else lam$P
  )
}

