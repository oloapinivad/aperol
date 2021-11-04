library("ggplot2")
library("cowplot")

# functions for the tuning
regions_cols <- c("darkred", "darkgreen", "darkblue")


# create a function to recreate the 3d multilinear fit
# use sweep and sum, there is no smart way
refit <- function(para, ml, square = "linear") {
  if (square == "quadratic") {
    paralist <- c(para, para^2)
    if (length(para) > 1) {
      interaction_term <- apply(combn(para, 2), 2, prod)
      paralist <- c(paralist, interaction_term)
    }
  } else if (square == "linear") {
    paralist <- para
  } else if (square == "mixed") {
    paralist <- c(para, para["GFRCRIT"]*para["GKDRAG"], para["GFRCRIT"]*para["GKWAKE"])
  } else if (square == "mixed2") {
    paralist <- c(para, para["GFRCRIT"]*para["GKDRAG"], para["GFRCRIT"]*para["GKWAKE"], para["ZTOFD"]*para["GKWAKE"]) 
  }
  out <- ml[1, , ] + apply(sweep(ml[-1, , , drop = F], 1, paralist, "*"), 2:3, sum)
  return(out)
}

# root mean square error computation with ww
rmse_core <- function(field, ww) {
  out <- sqrt(weighted.sum(
    field^2, ww
  ))
  return(out)
}

# declare rmse function from multilinear fit
rmse_level <- function(para, ml, lons, lats, ww, standard = F, square = "linear") {
  kk <- refit(para, ml, square)
  if (standard) {
    kk <- minmax_standardize(kk)
  }
  rmse <- rmse_core(kk[lons, lats], ww[lons, lats])
  return(rmse)
}


# main function that loops on levels
rmse_total <- function(para, para_exclude, region,
                       ww, multilinear, residuals,
                       residual_weights = F, standard = F,
                       square = "linear") {
  rmse <- 0
  for (season in seasons) {
    for (var in vars) {
      for (level in levels) {
        rr <- sel_region(region)
        ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]
        if (residual_weights) {
          ww <- ww / sse[, , which(var == vars), which(level == levels), which(season == seasons)]
        }
        paras <- c(para, para_exclude)
        rmse <- rmse + rmse_level(paras[order(names(paras))], ml, rr$lons, rr$lats, ww, standard = standard, square = square)
      }
    }
  }
  return(rmse)
}


# minmax standardize
minmax_standardize <- function(timeseries) {
  out <- (timeseries - min(timeseries, na.rm = T)) / (max(timeseries, na.rm = T) - min(timeseries, na.rm = T))
  return(out)
}


# normalize params
eval_norma <- function(params_list, centering = T) {
  mean_params <- colMeans(params_list)
  sd_params <- apply(params_list,2,sd)
  if (!centering) {
    mean_params <- mean_params * 0
    sd_params <- sd_params / sd_params
  }
  return(list(mean = mean_params, sd = sd_params, centering = centering))
}

# denormalize and normalize
denorma <- function(params, norma_object) {
  selected <- names(params)
  out <- params * norma_object$sd[selected] + norma_object$mean[selected]
  return(out)
}

norma <- function(params, norma_object) {
  selected <- names(params)
  out <- (params - norma_object$mean[selected]) / norma_object$sd[selected]
  return(out)
}








  
