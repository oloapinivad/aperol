# ensemble means
MILESDIR <- "/home/paolo/MiLES"
FILEDIR <- "/work/users/paolo/data"
TUNEDIR <- "/home/paolo/RScript/reforge/tune_gwd/aggressive_tuning"
source(file.path(MILESDIR, "script", "basis_functions.R"))
source(file.path(TUNEDIR, "..", "..", "reforge_help_functions.R"))
source(file.path(TUNEDIR, "aggressive_functions.R"))

# fundamental control flags
do_external <- T
do_loading <- T
do_compute <- T


#### ----- FLAGS TO BE CONTROLLED EXTERNALLY -----####
# variables on which you want to try the optimization
if (do_external) {
  args <- commandArgs(TRUE)
  names_args <- c(
    "vars", "nparams", "params_exclude", "levels",
    "sqorder", "seasons", "do_centering"
  )
  print("--------------------")
  print("READING COMMAND ARGS")
  for (k in seq_along(names_args)) {
    if (names_args[k] %in% c("levels", "nparams", "seasons")) {
      value <- eval(parse(text = args[k]))
    } else {
      value <- as.character(parse(text = args[k]))
    }

    assign(names_args[k], value)
    print(paste(names_args[k], paste(get(names_args[k]), collapse = "-")))
  }
  print("--------------------")
}

# do you want to weight with R-Squared the RMSE
do_residual_weights <- F

# flag for standardization
do_standard <- F

# regions
regions <- c("global", "northern_hemisphere")


if (!do_external) {

  # do you want to center the predictors
  do_centering <- T

  # nparams <- c("GKDRAG", "GFRCRIT")
  nparams <- c("ZTOFD", "GKWAKE", "GKDRAG", "GFRCRIT")

  # flag for 2nd order polynomial fit (false for linear fit)
  # sqorder <- "quadratic"
  sqorder <- "linear"
  #sqorder <- "mixed"
  #sqorder <- "mixed2"

  levels <- c(1000, 5000, 10000, 20000, 30000, 50000, 70000, 85000, 92500)
  levels <- c(5000, 25000, 85000)
  #levels <- c(85000)

  # seasons to be used together
  seasons <- "DJFM"
  # seasons <- c("DJFM", "JJAS")

  params_exclude <- c("GFRCRIT")
  # params_exclude <- "NULL"

  # variables
  vars <- c("ua")
}


##### --- END OF FLAG TO BE CONTROLLED EXTERNNALY -------###

# params included in the optimizations
params_include <- nparams[!nparams %in% params_exclude]

# deatils
project <- "REFORGE"
model <- "EC-Earth3"

# experiments used and ensemble available
expnames <- c("rfrg-ctrl-param", "rfrg-ctrl-sobol")

# define available ensembles from loadens2() function (add one for ctrl-param)
nensembles <- length(loadens2("EC-Earth3", expnames[2])) + 1

# use also default run to tune the model? deprecated flag
extra_default <- F

# convergene methods for optimizations
# conv_methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
# conv_methods <- c("Nelder-Mead")
conv_methods <- c("L-BFGS-B") # it allows for defining boundaries
if (length(nparams) == 1) {
  conv_methods <- "Brent"
}

# years
year1 <- 1999
year2 <- 2005

# checks
if (params_exclude == "NULL") {
  params_exclude <- NULL
}
do_centering <- as.logical(do_centering)


# folder destionation
kvars <- paste(paste(rep(vars, each = length(levels)), levels, sep = ""), collapse = "-")
kparams <- paste(nparams, collapse = "-")
kparams_exclude <- ifelse(is.null(params_exclude), "None", paste(params_exclude, collapse = "-"))
kresidual <- ifelse(do_residual_weights, ifelse(do_standard, "Residual+Standard", "Residual"), ifelse(do_standard, "Standard", "Default"))
kseason <- paste(seasons, collapse = "-")
kcentering <- ifelse(do_centering, "Centering", "Original")

# sort params
nparams <- sort(nparams)
params_include <- sort(params_include)
# params_exclude <- sort(params_exclude)
print(sqorder)

# properties for quadratic or linear fit
if (sqorder == "quadratic") {
  params_dim <- 1 + length(nparams) * 2 + length(nparams) * (length(nparams) - 1) / 2
  # squared as combination of I(PARAM^2) terms and linear terms (use .^2 to get interaction terms)
  regre_formula <- paste(
    "y ~ .^2", " + ",
    paste0(rep("I(", length(nparams)), nparams, rep("^2)", length(nparams)), collapse = " + ")
  )
} else if (sqorder == "linear") {
  params_dim <- length(nparams) + 1
  regre_formula <- paste("y ~", paste0(nparams, collapse = " + "))
} else if (sqorder == "mixed") {
  params_dim <- length(nparams) + 3
  regre_formula <- paste("y ~", paste0(nparams, collapse = " + "), paste0("+ GKDRAG:GFRCRIT + GKWAKE:GFRCRIT"))
  # regre_formula <- paste("y ~", paste0("GKDRAG:GFRCRIT + GKWAKE:GFRCRIT + "), paste0(sample(nparams,4), collapse = " + "))
} else if (sqorder == "mixed2") {
  params_dim <- length(nparams) + 4
  regre_formula <- paste("y ~", paste0(nparams, collapse = " + "), paste0("+ GKDRAG:GFRCRIT + GKWAKE:GFRCRIT + GKWAKE:ZTOFD"))
}



# number of exp to loop on
if (do_external) {
  nelements <- (params_dim + 2):nensembles
  nloop <- 19
} else {
  nelements <- nensembles
  nloop <- 1
}

# figure dire
# create folder
FIGDIR <- file.path(
  "/work/users/paolo/figures/REFORGE/tune_gwd/new_fingerprint/",
  paste0(year1, "-", year2), kvars, kparams, kparams_exclude, kcentering, sqorder, kresidual, kseason
)
dir.create(FIGDIR, recursive = T)

# pvalue for the multiinear regression
PVALUE <- 0.05


# file loading
if (do_loading) {

  # dataframe for expeirments RMSE
  expdf <- data.frame()
  for (season in seasons) {
    for (var in vars) {
      for (level in levels) {
        print(level)

        # loading of ERA5
        filename <- Sys.glob(file.path(FILEDIR, "ERA5", field_freq(var), var, "*.nc"))
        print(filename)
        field <- ncdf.opener.universal(filename,
          namevar = var, tmonths = season2timeseason(season),
          tlev = loadlevel(level), verbose = F,
          tyears = 1979:2019
        )
        era_field <- rowMeans(field$field, dim = 2)
        assign(paste0(var, level, "_ERA5_", season), era_field)

        # produce ww
        area_ww <- area_weights(field$lon, field$lat)

        # declare arrays for variablers and parameters
        new_array <- array(NA, dim = c(length(field$lon), length(field$lat), nensembles))
        params_complete <- array(NA, dim = c(nensembles, length(nparams)), dimnames = list(NULL, nparams))

        # loading model data
        count <- 0
        for (expname in expnames) {
          for (ensemble in loadens2(model, expname)) {
            filename <- Sys.glob(file.path(
              FILEDIR, "REFORGE", model, expname, ensemble,
              field_freq(var), var, "*.nc"
            ))

            # load files
            if (length(filename) > 0) {
              print(filename)
              if (expname == "rfrg-ctrl-param") {
                years <- 2000:2029
              } else {
                years <- year1:year2
              }
              field <- ncdf.opener.universal(filename,
                namevar = var, tmonths = season2timeseason(season),
                tlev = loadlevel(level), tyears = years
              )

              # assign loaded file and save its bias with respect to era5
              mean_field <- rowMeans(field$field, dim = 2)
              ens_info <- unlist(ensname(expname, ensemble)$params)

              # handle the extra ensembles, a weird shytti code
              count <- count + 1
              for (nparam in nparams) {
                params_complete[count, nparam] <- ens_info[nparam]
              }
              # params_complete[count, ] <- ens_info[names(ens_info) %in% nparams]
              new_array[, , count] <- mean_field - era_field

              # estimate immediately its own rmse against ERA5 in the required regions
              # put everything into a dataframe
              for (region in regions) {
                rr <- sel_region(region)
                rmse <- rmse_core(
                  new_array[rr$lons, rr$lats, count],
                  area_ww[rr$lons, rr$lats]
                )
                # assign value of each parameters for each ensemble
                for (nparam in nparams) {

                  # get values from the params arrary
                  v <- params_complete[count, nparam]

                  # compute wind speed
                  wind_speed <- mean(abs(mean_field[rr$lons, rr$lats]))

                  # build up the dataframe
                  element <- data.frame(
                    exp = expname, ens = ensemble, variable = var, lev = level,
                    parameter = nparam, paramvalue = v[[1]],
                    reg = region, rmse = rmse, wind_speed = wind_speed, season = season
                  )
                  expdf <- rbind(expdf, element, stringsAsFactors = F)
                }
              }
            }
            # assign the bias array
            assign(
              paste0(
                var, level, "_", model, "_new_array_", season
              ),
              new_array
            )

            assign(
              paste0(
                var, level, "_", model, "_",
                expname, "_", ensemble, "_", season
              ),
              mean_field - era_field
            )
            # print(mean(mean_field - era_field))
          }
        }
      }
    }
  }
  # axes
  lon <- field$lon
  lat <- field$lat

  # save the params values
  params_original <- params_complete
}

# core business of computation
if (do_compute) {
  swapping <- evolution <- data.frame()

  # normalization, applied only if do_centering=T
  nor <- eval_norma(params_original, centering = do_centering)
  params_complete <- sweep(sweep(params_original, 2, nor$mean), 2, nor$sd, "/")

  # if not include the control experiment
  if (!extra_default) {
    nelements <- nelements - 1
  }

  # loop to sample the solution on the ensembles and on the swapping possibilities
  for (nelement in nelements) {
    print(nelement)
    print(nelements)
    for (loop in 1:nloop) {
      if (nelement == (max(nelements) - 1)) {
        samples <- (1:max(nelements))[-loop]
      } else if (nelement == max(nelements)) {
        samples <- 1:nelement
        if (loop >= 2) {
          next
        }
      } else {
        if (loop == 1) {
          samples <- 1:nelement
        } else {
          samples <- sample(max(nelements), nelement)
        }
      }
      print(samples)

      revolution <- NULL
      print(paste("-------> Loop on only", nelement, "ensembles with nloop...", loop))


      # multilinear regression
      print("Multinear regression...")

      pvalue <- multilinear <- explained <- array(NA, dim = c(params_dim, length(lon), length(lat), length(vars), length(levels), length(seasons)))
      sse <- rsquared <- array(NA, dim = c(length(lon), length(lat), length(vars), length(levels), length(seasons)))
      for (season in seasons) {
        for (var in vars) {
          # for (level in field_levels(var)) {
          for (level in levels) {
            print(paste(season, var, level))
            complete_array <- get(paste0(var, level, "_", model, "_new_array_", season))

            # remove ctrl experiment
            if (!extra_default) {
              complete_array <- complete_array[, , -1]
              complete_params <- params_complete[-1, , drop = F]
            }

            # load only sampled experiments
            new_array <- complete_array[, , samples]
            params <- complete_params[samples, , drop = F]

            # multilinear regression, for the moment avoid standardization of values - apply standardization to RMSE calculation
            # new_array <- standardize(new_array)

            # produce a formula as a string linear combination of the params: use data.frame to simplify the names
            params <- as.data.frame(params)

            # define the massive linear fit as ndimensional list
            print("Big regression...")
            beast <- apply(new_array, 1:2, function(y) {
              lm(as.formula(regre_formula), data = params)
            })

            # beast2 <- apply(new_array, 1:2, function(y) {
            #  glm(as.formula(regre_formula), data = params)
            # })

            # save dimension names
            coeff_names <- dimnames(sapply(beast, "[[", "coefficients"))[[1]]

            # create the summary ndimensional list and then extract rsquared
            print("Rsquared...")
            summary_beast <- lapply(beast, summary)
            rq <- sapply(summary_beast, "[[", "adj.r.squared")
            rsquared[, , which(var == vars), which(level == levels), which(season == seasons)] <- rq
            revolution <- c(revolution, mean(rq))

            # apply anova on the fit and extract variance explained and the sse
            print("Var_explained...")
            anova_beast <- lapply(beast, anova)
            tmp_explained <- array(sapply(anova_beast, "[[", "Sum Sq"), dim = c(params_dim, length(lon), length(lat)))
            sse[, , which(var == vars), which(level == levels), which(season == seasons)] <- tmp_explained[params_dim, , ]
            explained[, , , which(var == vars), which(level == levels), which(season == seasons)] <-
              sweep(tmp_explained, 2:3, apply(tmp_explained, 2:3, sum), "/")
            # issue in dimension naming: residuals are the last in the anova!
            dimnames(explained)[[1]] <- c(coeff_names[2:length(coeff_names)], coeff_names[1])

            # extract coefficients
            print("Coefficients...")
            multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)] <-
              sapply(beast, "[[", "coefficients")
            dimnames(multilinear)[[1]] <- coeff_names

            # verify the rebuild quality of the other elements
            print("Checking rebuilding experiments...")
            excluded_elements <- which(!1:max(nelements) %in% samples)
            ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]
            for (excluded_element in excluded_elements) {
              print(paste("Excluded element number:", excluded_element))
              parameter_set <- unlist(complete_params[excluded_element, ])
              rebuild_field <- refit(parameter_set, ml, square = sqorder)
              original_field <- complete_array[, , excluded_element]
              delta_field <- rebuild_field - original_field
              for (region in regions) {
                rr <- sel_region(region)
                val <- rmse_core(delta_field[rr$lons, rr$lats], area_ww[rr$lons, rr$lats])
                subelement <- data.frame(
                  region = region, nens = excluded_element,
                  rmse = val, var = var, level = level, season = season,
                  loop = loop, nelement = nelement
                )
                swapping <- rbind(swapping, subelement, stringsAsFactors = FALSE)
              }
            }

            assign(paste0(var, level, "_", model, "_beast_", season), beast)
            assign(paste0(var, level, "_", model, "_summary_beast_", season), summary_beast)
            assign(paste0(var, level, "_", model, "_anova_beast_", season), anova_beast)
          }
        }
      }

      # define the optimization functions!
      print("Optimizations...")

      # optimization on the two main regions

      # default for excluded params and mean values for included ones to adjust boundaries of convergence
      if (is.null(params_exclude)) {
        param_default_exclude <- NULL
      } else {
        param_default_exclude <- norma(sapply(params_exclude, function(x) params_properties(x)$pdef), nor)
      }

      # extract default params with sapply
      param_default_include <- norma(sapply(params_include, function(x) params_properties(x)$pdef), nor)


      # loop on convergence methods and number of tentatives
      for (region in regions) {
        for (conv_method in conv_methods) {
          print(conv_method)
          if (conv_method %in% c("L-BFGS-B", "Brent")) {
            param_lower <- norma(sapply(params_include, function(x) params_properties(x)$pmin), nor)
            param_upper <- norma(sapply(params_include, function(x) params_properties(x)$pmax), nor)
          } else {
            param_lower <- -Inf
            param_upper <- Inf
          }

          # call the real optimizations with all the flags
          best <- optim(param_default_include, rmse_total,
            para_exclude = param_default_exclude,
            region = region,
            ww = area_ww,
            multilinear = multilinear,
            residual_weights = do_residual_weights,
            residuals = sse,
            standard = do_standard,
            square = sqorder,
            method = conv_method,
            lower = param_lower,
            upper = param_upper
          )
          # print(best$value)
          # print(best$par)
        }

        # reintegrated removed par and keep correct order
        best$par <- c(best$par, param_default_exclude)
        best$par <- best$par[order(names(best$par))]
        print(best$par)

        if (length(nparams) == 1) {
          names(best$par) <- nparams
        }

        # save evolution in dataframe
        element <- data.frame(
          region = region, nens = nelement,
          rmse = best$value, rsquared = mean(revolution), loop = loop
        )
        element <- cbind(element, as.data.frame(t( denorma(best$par, nor))))
        print(element)
        evolution <- rbind(evolution, element, stringsAsFactors = F)
        assign(paste0("best_", region), best)
      }
    }
  }

  print("Complete computation...")
  for (season in seasons) {
    for (var in vars) {
      # for (level in field_levels(var)) {
      for (level in levels) {
        print(paste(season, var, level))

        beast <- get(paste0(var, level, "_", model, "_beast_", season))
        summary_beast <- get(paste0(var, level, "_", model, "_summary_beast_", season))

        print("P-values...")
        pvalue[, , , which(var == vars), which(level == levels), which(season == seasons)] <-
          sapply(lapply(summary_beast, "[[", "coefficients"), "[", , 4)
        dimnames(pvalue)[[1]] <- coeff_names
      }
    }
  }

  # test reconstruction with the multilinear regressions making use of
  # set of permutations of parametrs and then averaging the results
  # as a function of the parameter of interest: check for RMSE!
  # print("Permutations...")
  # newdf <- data.frame()
  # for (season in seasons) {
  #  for (var in vars) {
  #    # for (level in field_levels(var)) {
  #    for (level in levels) {
  #      print(paste(var, level))
  #      ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]
  #
  #        # sapply and lapply to create a list to be expanded from nparams
  #        param_to_permute <- sapply(nparams, function(x) params_properties(x)$pseq)
  #        param_to_list <- lapply(seq_len(ncol(param_to_permute)), function(i) param_to_permute[, i])
  #        names(param_to_list) <- nparams
  #
  #        # expand grid to create the permutations
  #        permutations <- expand.grid(param_to_list)
  #        subperm <- sample(length(permutations[, 1]), min(length(permutations[, 1]), 100))
  #        # subperm <- seq_along(permutations[, 1])
  #        for (p in subperm) {

  # create the field for the parameter set
  #          parameter_set <- unlist(permutations[p, ])
  #          rebuild_field <- refit(parameter_set, ml, square = sqorder)
  #
  #          # estimate rmse
  #          for (region in regions) {
  #            rr <- sel_region(region)
  #            val <- rmse_core(rebuild_field[rr$lons, rr$lats], area_ww[rr$lons, rr$lats])
  #
  #            # computed full field
  #            abs_field <- rebuild_field + get(paste0(var, level, "_ERA5_", season))
  #
  # extract northern hemisphere wind, for Atlantic and Global
  #            wind_speed <- mean(abs_field[rr$lons, whicher(field$lat, 20):whicher(field$lat, 80)])
  #
  #            # create storage dataframe
  #            for (nparam in nparams) {
  #              element <- data.frame(
  #                parameter = nparam, variable = var, lev = level,
  #                paramvalue = parameter_set[which(nparam == nparams)][[1]],
  #                reg = region, rmse = val, wind_speed = wind_speed, season = season
  #              )
  #              newdf <- rbind(newdf, element, stringsAsFactors = F)
  #            }
  #          }
  #        }
  #      }
  #    }
  #  }

  # evaluate the capacity of rebuild current experiments
  # print("Rebuilding...")
  # rebuildf <- data.frame()
  # for (season in seasons) {
  #    for (var in vars) {
  #      # for (level in field_levels(var)) {
  #      for (level in levels) {
  #        print(paste(season, var, level))
  #        ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]
  #
  #        for (p in seq_along(params_complete[, 1])) {
  #          # create the field for the parameter set
  #          parameter_set <- unlist(params_complete[p, ])
  #          rebuild_field <- refit(parameter_set, ml, square = sqorder)
  #
  #          # estimate rmse
  #          for (region in regions) {
  #            rr <- sel_region(region)
  #            val <- rmse_core(rebuild_field[rr$lons, rr$lats], area_ww[rr$lons, rr$lats])

  #            # computed full field
  #            abs_field <- rebuild_field + get(paste0(var, level, "_ERA5_", season))
  #
  #            # extract northern hemisphere wind, for Atlantic and Global
  #            wind_speed <- mean(abs_field[rr$lons, whicher(field$lat, 20):whicher(field$lat, 80)])
  #
  #            # create storage dataframe
  #            for (nparam in nparams) {
  #              element <- data.frame(
  #                rebuild_rmse = val, rebuild_wind_speed = wind_speed
  #              )
  #              rebuildf <- rbind(rebuildf, element, stringsAsFactors = F)
  #            }
  #          }
  #        }
  #      }
  #    }
  #  }
  #
  #  expdf <- cbind(expdf, rebuildf)


  # summarize the results of RMSE
  # print("Summary...")
  # summarydf <- summary.dataframe(newdf[, !names(newdf) %in% "wind_speed"],
  #  measurevar = "rmse", groupvars = c("variable", "lev", "reg", "paramvalue", "parameter", "season")
  # )
  # summarywind <- summary.dataframe(newdf[, !names(newdf) %in% "rmse"],
  #  measurevar = "wind_speed", groupvars = c("variable", "lev", "reg", "paramvalue", "parameter", "season")
  # )
}


# print pvalue coefficients single parameter plots
print("Plotting...")
onedims <- c("Rsquared", "SSE")
for (onedim in onedims) {
  print(onedim)
  name <- paste(FIGDIR, "/", onedim, "_",
    year1, "-", year2, ".pdf",
    sep = ""
  )
  pdf(
    file = name, width = 10 * length(levels), height = 6 * length(vars) * length(seasons),
    onefile = T, bg = "white", family = "Helvetica"
  )
  panels <- c(length(seasons) * length(vars), length(levels))
  par(c(plotpar, list(mfrow = panels)))

  for (season in seasons) {
    count <- 0
    for (var in vars) {
      for (level in levels) {
        fp <- field_properties(var, level)
        if (onedim == "Rsquared") {
          fp$color_diff <- palette.gold
          fp$lev_diff <- seq(0, 1, 0.1)
          plot_diff <- rsquared[, , which(var == vars), which(level == levels), which(season == seasons)]
        } else if (onedim == "SSE") {
          fp$color_diff <- palette.spct
          fp$lev_diff <- seq(0, 100, 5)
          plot_diff <- sse[, , which(var == vars), which(level == levels), which(season == seasons)]
        }

        print(paste(season, var, level))
        count <- count + 1
        plot_title <- paste(onedim, var, level, season)
        plot_full <- get(paste0(var, level, "_ERA5_", season))

        im1 <- plot.prepare(lon, lat, plot_diff,
          proj = map_projection, lat_lim = lat_lim
        )
        im2 <- plot.prepare(lon, lat, plot_full * fp$factor,
          proj = map_projection, lat_lim = lat_lim
        )
        filled.contour3(im1$x, im1$y, im1$z,
          xlab = im1$xlab, ylab = im1$ylab, main = plot_title,
          levels = fp$lev_diff, color.palette = fp$color_diff, asp = im1$asp, extend = T,
          xlim = im1$xlim, ylim = im1$ylim, axes = im1$axes
        )
        for (region in regions) {
          rr <- sel_region(region)
          expvar <- mean(plot_diff[rr$lons, rr$lats])
          text(100, 90 - which(region == regions) * 10, paste0(rr$sname, onedim, ":", round(expvar, 2)),
            cex = 2.5, col = regions_cols[which(region == regions)]
          )
        }

        contour(im2$x, im2$y, im2$z,
          levels = fp$lev_field,
          add = T, lwd = 2, drawlabels = T, labcex = clab
        )
        contour(im2$x, im2$y, im2$z,
          levels = -fp$lev_field,
          add = T, lwd = 2, drawlabels = T, lty = 2, labcex = clab
        )
        proj.addland(lon, lat, proj = map_projection)
        mtext(lettering[count], line = 1, adj = 0, cex = cex.letter)
        if (count == length(vars) * length(seasons) * length(levels)) {
          image.scale3(volcano,
            levels = fp$lev_diff, color.palette = fp$color_diff,
            colorbar.label = fp$legend_unit,
            cex.colorbar = imgscl_colorbar, cex.label = imgscl_label,
            colorbar.width = 1.5, line.colorbar = 0,
            line.label = fp$legend_distance
          )
        }
      }
    }
  }
  dev.off()
}

# print regression coefficients and strippling with pvalue
twodims <- c("Regression_coefficients", "Explained_variance")

for (twodim in twodims) {
  print(twodim)
  for (season in seasons) {
    for (var in vars) {
      kp <- c(length(levels), params_dim - 1)
      name <- paste(FIGDIR, "/", twodim, "_", var, "_",
        year1, "-", year2, "_", season, ".pdf",
        sep = ""
      )
      pdf(
        file = name, width = 8 * kp[2], height = 5 * kp[1],
        onefile = T, bg = "white", family = "Helvetica"
      )
      par(c(plotpar, list(mfrow = kp)))
      for (level in levels) {
        print(paste(season, var, level))
        fp <- field_properties(var, level)

        if (twodim == "Regression_coefficients") {
          ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]
          pv <- pvalue[, , , which(var == vars), which(level == levels), which(season == seasons)]
          fp$color_diff <- palette.rdbl
        } else if (twodim == "Explained_variance") {
          ml <- explained[, , , which(var == vars), which(level == levels), which(season == seasons)]
          fp$color_diff <- palette.bupu
          lev_diff <- seq(0, 1, 0.1)
          plegend <- paste("Explained variance")
        }

        count <- 0
        for (coeff in coeff_names[-1]) {
          if (twodim == "Regression_coefficients") {
            if (do_centering) {
              lev_diff <- fp$lev_diff / 4
              extralegend <- paste("per normalized", coeff)

            } else { 
              lev_diff <- fp$lev_diff / mean(range(params_properties(coeff)$pseq))
              extralegend <- paste( "per", coeff, "unit")
            }
            plot_pvalue <- pv[coeff, , ]
            plot_pvalue <- ifelse(plot_pvalue <= PVALUE, 1, NA)
            plegend <- paste(fp$legend_unit, extralegend)
          }

          plot_diff <- ml[coeff, , ]
          count <- count + 1
          plot_full <- get(paste0(var, level, "_ERA5_", season))

          plot_title <- paste(coeff, var, level, season)

          im1 <- plot.prepare(lon, lat, plot_diff * fp$factor,
            proj = map_projection, lat_lim = lat_lim
          )
          im2 <- plot.prepare(lon, lat, plot_full * fp$factor,
            proj = map_projection, lat_lim = lat_lim
          )
          filled.contour3(im1$x, im1$y, im1$z,
            xlab = im1$xlab, ylab = im1$ylab, main = plot_title,
            levels = lev_diff, color.palette = fp$color_diff, asp = im1$asp, extend = T,
            xlim = im1$xlim, ylim = im1$ylim, axes = im1$axes
          )
          contour(im2$x, im2$y, im2$z,
            levels = fp$lev_field,
            add = T, lwd = 2, drawlabels = T, labcex = clab
          )
          contour(im2$x, im2$y, im2$z,
            levels = -fp$lev_field,
            add = T, lwd = 2, drawlabels = T, lty = 2, labcex = clab
          )
          if (twodim == "Regression_coefficients") {
            im3 <- plot.prepare(lon, lat, plot_pvalue,
              proj = map_projection, lat_lim = lat_lim
            )
            xx <- expand.grid(x = im3$x, y = im3$y)
            points(xx$x[im3$z > 0.5], xx$y[im3$z > 0.5], pch = 20, cex = 0.5)
          } else if (twodim == "Explained_variance") {
            for (region in regions) {
              rr <- sel_region(region)
              expvar <- mean(plot_diff[rr$lons, rr$lats])
              text(100, 90 - which(region == regions) * 10, paste0(rr$sname, "ExpVar:", round(expvar, 2)),
                cex = 2.5, col = regions_cols[which(region == regions)]
              )
            }
          }

          proj.addland(lon, lat, proj = map_projection)
          mtext(lettering[count], line = 1, adj = 0, cex = cex.letter)
          if (twodim == "Explained_variance" & count < kp[2] * kp[1]) {
            next
          }
          image.scale3(volcano,
            levels = lev_diff, color.palette = fp$color_diff,
            colorbar.label = plegend,
            cex.colorbar = imgscl_colorbar, cex.label = imgscl_label,
            colorbar.width = 1.5, line.colorbar = 1,
            line.label = fp$legend_distance
          )
        }
      }
      dev.off()
    }
  }
}


# block to evaluate the quality of the RMSE convergence compared to the other experiments
print("RMSE convergence...")
check <- summary.dataframe(subset(expdf, parameter == nparams[1]),
  measurevar = "rmse", groupvar = c("exp", "reg", "ens")
)[, c("exp", "reg", "rmse")]
check$rmse <- check$rmse * length(levels) * length(vars)
final <- subset(evolution, nens == (nensembles - 1), select = c("region", "rmse", "nens"))
colnames(final) <- c("reg", "rmse", "exp")
convergencedf <- rbind(final, check)
theplot <- ggplot(convergencedf, aes(x = reg, y = rmse, fill = exp, col = exp)) +
  geom_boxplot() +
  labs(
    title = paste("RMSE from different experiments"),
    x = paste("Region"), y = "RMSE (m/s)", col = "Experiments", fill = "Experiments"
  ) +
  stat_summary(mapping = aes(label = sprintf("%0.3f", ..y..)), fun = mean, geom = "text", vjust = -0.25, size = 6) +
  coord_cartesian(ylim = c(0, 3) * length(levels) * length(seasons)) +
  theme_light(base_size = 18)

name <- file.path(FIGDIR, paste(
  "RMSE_compared_to_exps.pdf",
  sep = ""
))
save_plot(name, theplot, base_width = 6, base_height = 6)


# error in rebuilding missing experiments
if (length(swapping)>0) {
if (length(swapping[,1]) > 8) {
  print("RMSE rebuilding quality...")
  theplot <- list()
  count <- 0
  for (s in seasons) {
    for (v in vars) {
      for (l in levels) {
        count <- count + 1
        mydf <- subset(swapping, season == s & level == l & var == v)
        mydf$grouper <- interaction(mydf$nelement, mydf$region)
        theplot[[count]] <- ggplot(mydf, aes(x = nelement, y = rmse, fill = region)) +
          geom_boxplot(aes(group = grouper)) +
          labs(
            title = paste("Rebuild RMSE", v, l, s),
            x = paste("# members"), y = "RMSE", col = "Region"
          ) +
          stat_summary(mapping = aes(col = region), fun = mean, geom = "line", size = 2) +
          stat_summary(mapping = aes(label = sprintf("%0.2f", ..y..)), fun = mean, geom = "text", vjust = -0.25, size = 6) +
          coord_cartesian(xlim = c(6, 20), ylim = c(0, 2)) +
          theme_light(base_size = 18)
        if (count == length(levels) * length(vars) * length(seasons)) {
          plot_legend <- get_legend(theplot[[count]])
        }
        theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
      }
    }
  }
  prow <- plot_grid(
    plotlist = theplot, align = "vh",
    labels = lettering[1:count], hjust = 0,
    ncol = length(levels), label_size = 24
  )
  name <- file.path(FIGDIR, paste(
    "RMSE_rebuilding_quality.pdf",
    sep = ""
  ))
  prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 2), nrow = 1, ncol = 2)
  save_plot(name, prow2, base_width = 6 * length(levels), base_height = 6 * length(seasons))
}
}

# plot for RMSE sensitivity to each parameter
# show the mean value from the permutation excercise, ribbons as standard deviation
# but also the value from the real experiments to see if they fit
print("Sensitivity to parameters...")
for (kind in c("rmse", "wind_speed")) {
  for (region in regions) {
    theplot <- list()
    count <- 0
    for (seas in seasons) {
      for (nparam in nparams) {
        for (var in vars) {
          fp <- field_properties(var, level)
          count <- count + 1
          if (kind == "rmse") {
            LL <- c(0, 2) * mean(abs(fp$lev_diff)) * fp$factor
          } else if (kind == "wind_speed") {
            LL <- c(0, 8) * fp$factor * mean(abs(fp$lev_diff))
          }
          subexpdf <- subset(subset(expdf, season == seas), var == variable & parameter == nparam & reg == region)
          theplot[[count]] <- ggplot(subexpdf, mapping = aes_string(
            x = "paramvalue", y = kind,
            col = "factor(lev)"
          )) +
            geom_point() +
            geom_smooth(method = "lm", formula = y ~ x) +
            # scale_y_log10() +
            labs(
              title = paste(region, var, nparam), x = paste(nparam),
              y = kind, fill = "Level (Pa)", col = "Level (Pa)"
            ) +
            coord_cartesian(ylim = LL) +
            theme_light(base_size = 18)

          if (count == length(nparams) * length(vars)) {
            plot_legend <- get_legend(theplot[[count]])
          }
          theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
        }
      }
    }
    lettering <- paste0("(", letters, ")")
    name <- file.path(FIGDIR, paste(
      kind, "_", region, "_sensitivity_to_params.pdf",
      sep = ""
    ))
    prow <- plot_grid(
      plotlist = theplot, align = "vh",
      labels = lettering[1:count], hjust = 0,
      ncol = length(nparams), label_size = 24
    )
    prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 1.1), nrow = 1, ncol = 2)
    save_plot(name, prow2, base_width = 6 * length(nparams), base_height = 6 * length(vars))
  }
}

# plot for RMSE sensitivity to each parameter
# compare rebuilt with original exps
# kinds <- "rmse"
# for (season in seasons) {
#  for (kind in kinds) {
#    print(kind)
#    for (region in regions) {
#      theplot <- list()
#      count <- 0
#      for (nparam in nparams[1]) {
#        for (var in vars) {
#          fp <- field_properties(var, level)
#          if (kind == "rmse") {
#            LL <- c(0, 5) * mean(abs(fp$lev_diff)) / fp$factor
#          } else if (kind == "wind_speed") {
#            LL <- c(0, 20) * mean(abs(fp$lev_diff)) / fp$factor
#          }
#
#          count <- count + 1
#          # for (level in field_levels(var)) {
#          subexpdf <- subset(expdf, var == variable & parameter == nparam & reg == region)
#          theplot[[count]] <- ggplot(subexpdf) +
#            geom_point(mapping = aes_string(
#              x = paste0("rebuild_", kind), y = kind,
#              col = "factor(lev)", shape = "exp", size = 2
#            )) +
#            # scale_y_log10() +
#            labs(title = paste(region, var), x = paste("Estimated", kind), y = kind, shape = "Expname", col = "Level (Pa)") +
#            coord_cartesian(xlim = LL, ylim = LL) +
#            theme_light(base_size = 18)
#
#          if (count == length(nparams) * length(vars)) {
#            plot_legend <- get_legend(theplot[[count]])
#          }
#          theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
#        }
#      }
#      lettering <- paste0("(", letters, ")")
#      name <- paste(FIGDIR, "/",
#        kind, "_comparisons_to_original_params_", region, "_", season, ".pdf",
#        sep = ""
#      )
#      prow <- plot_grid(
#        plotlist = theplot, align = "vh",
#        labels = lettering[1:count], hjust = 0,
#        ncol = length(vars), label_size = 24
#      )
#      #prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 1.1), nrow = 1, ncol = 2)
#      save_plot(name, prow, base_width = 10 * length(vars) + 2, base_height = 10)
#    }
#  }
# }

# Plot the multilinear regression for some specific parameters set (as default)
# and the optimized one, which might be the results from the final convergence of parameters
for (season in seasons) {
  for (var in vars) {
    kp <- c(length(levels), 4)
    name <- paste(FIGDIR, "/EvaluateControl_",
      year1, "-", year2, "_", season, ".pdf",
      sep = ""
    )
    pdf(
      file = name, width = 10 * kp[2], height = 6 * kp[1],
      onefile = T, bg = "white", family = "Helvetica"
    )
    panels <- kp
    par(c(plotpar, list(mfrow = panels, cex.main = 2)))
    # plot properties

    for (level in levels) {
      ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]

      # add original plot for rfrg-ctrl-param
      plot_full <- get(paste0(var, level, "_ERA5_", season))

      # these defines the sets that are plotted, completely arbitrary
      sets <- c("rfrg-ctrl-param", "base", "delta", "skip")
      count <- 0
      for (set in sets) {
        fp <- field_properties(var, level)
        if (set == "rfrg-ctrl-param") {
          plot_diff <- get(paste0(
            var, level, "_", model, "_", set, "_",
            loadens2(model, set), "_", season
          ))
          plot_title <- set
        } else if (set == "skip") {
          if (level == levels[1]) {
            CX <- 3
            plot(0, 0, type = "n", axes = FALSE, ann = FALSE, xlim = c(0, 0.2), ylim = c(0, 0.2))
            text(0, 0.15, paste0("Vars: ", kvars), pos = 4, cex = CX)
            text(0, 0.125, paste0("Nensembles: ", nensembles - 1, " and Nyears:", length(year1:year2)), pos = 4, cex = CX)
            text(0, 0.1, paste0("knobs: ", kparams), pos = 4, cex = CX)
            text(0, 0.075, paste0("Fit type: ", sqorder), pos = 4, cex = CX)
            text(0, 0.05, paste0("Residual weighting:", do_residual_weights), pos = 4, cex = CX)
            text(0, 0.025, paste0("Standardize RMSE:", do_standard), pos = 4, cex = CX)
            text(0, 0, paste0("Centering:", do_centering), pos = 4, cex = CX)
          } else {
            plot.new()
          }
          next
        } else if (set == "base") {
          para <- params_complete[1, ]
          plot_diff <- refit(para, ml, square = sqorder)
          plot_title <- paste(set, "\n", paste0(nparams, ": ", round(denorma(para, nor), 2), collapse = " "))
        } else if (set == "delta") {
          plot_diff <- refit(para, ml, square = sqorder) - get(paste0(
            var, level, "_", model, "_", "rfrg-ctrl-param", "_",
            loadens2(model, set), "_", season
          ))
          plot_title <- set
          fp$color_diff <- palette.rdbl
        }

        # to recompose the linear fit
        count <- count + 1

        im1 <- plot.prepare(lon, lat, plot_diff * fp$factor,
          proj = map_projection, lat_lim = lat_lim
        )
        im2 <- plot.prepare(lon, lat, plot_full * fp$factor,
          proj = map_projection, lat_lim = lat_lim
        )
        filled.contour3(im1$x, im1$y, im1$z,
          xlab = im1$xlab, ylab = im1$ylab, main = plot_title,
          levels = fp$lev_diff, color.palette = fp$color_diff, asp = im1$asp, extend = T,
          xlim = im1$xlim, ylim = im1$ylim, axes = im1$axes
        )
        contour(im2$x, im2$y, im2$z,
          levels = fp$lev_field,
          add = T, lwd = 2, drawlabels = T, labcex = clab
        )
        contour(im2$x, im2$y, im2$z,
          levels = -fp$lev_field,
          add = T, lwd = 2, drawlabels = T, lty = 2, labcex = clab
        )
        kol <- c("darkred", "darkgreen", "darkblue")
        for (region in regions) {
          rr <- sel_region(region)
          # bias <- weighted.mean(plot_diff, ww)
          rmse <- sqrt(weighted.sum(plot_diff[rr$lons, rr$lats]^2, area_ww[rr$lons, rr$lats]))
          text(120, 90 - which(region == regions) * 10, paste0(rr$sname, "RMSE:", round(rmse, 2)),
            cex = 3, col = kol[which(region == regions)]
          )
        }

        # text(120, 70, paste0("RMSE:", round(rmse, 2)), cex = 3, col = "darkred")
        # text(120, 80, paste0("NA RMSE:", round(rmse_na, 2)), cex = 3, col = "darkgreen")
        proj.addland(lon, lat, proj = map_projection)
        mtext(lettering[count], line = 1, adj = -0.05, cex = cex.letter)
        image.scale3(volcano,
          levels = fp$lev_diff, color.palette = fp$color_diff,
          colorbar.label = fp$legend_unit,
          cex.colorbar = imgscl_colorbar, cex.label = imgscl_label,
          colorbar.width = 1.5, line.colorbar = 0,
          line.label = fp$legend_distance
        )
      }
    }
    dev.off()
  }
}




# Plot the multilinear regression for some specific parameters set (as default)
# and the optimized one, which might be the results from the final convergence of parameters
print("Multilinear regression...")
sets <- c("rfrg-ctrl-param", "base", paste0("best_", regions), "skip")
for (season in seasons) {
  for (var in vars) {
    fp <- field_properties(var, level)
    kp <- c(length(levels), length(sets))
    name <- paste(FIGDIR, "/", var, "_MultiLinearRegression_",
      year1, "-", year2, "_", season, ".pdf",
      sep = ""
    )
    pdf(
      file = name, width = 8 * kp[2], height = 5 * kp[1],
      onefile = T, bg = "white", family = "Helvetica"
    )
    panels <- kp
    par(c(plotpar, list(mfrow = panels, cex.main = 2)))
    for (level in levels) {
      print(paste(season, var, level))
      # plot properties

      ml <- multilinear[, , , which(var == vars), which(level == levels), which(season == seasons)]

      # add original plot for rfrg-ctrl-param
      plot_full <- get(paste0(var, level, "_ERA5_", season))
      #
      # these defines the sets that are plotted, completely arbitrary
      count <- 0
      for (set in sets) {
        if (set == "rfrg-ctrl-param") {
          plot_diff <- get(paste0(
            var, level, "_", model, "_", set, "_",
            loadens2(model, set), "_", season
          ))
          plot_title <- set
        } else if (set == "skip") {
          plot(0, 0, type = "n", axes = FALSE, ann = FALSE, xlim = c(0, 0.2), ylim = c(0, 0.2))
          if (level == levels[1]) {
            text(0, 0.15, paste0("Vars: ", kvars), pos = 4, cex = CX)
            text(0, 0.125, paste0("Nensembles: ", nensembles - 1, " and Nyears:", length(year1:year2)), pos = 4, cex = CX)
            text(0, 0.1, paste0("knobs: ", kparams), pos = 4, cex = CX)
            text(0, 0.075, paste0("Fit type: ", sqorder), pos = 4, cex = CX)
            text(0, 0.05, paste0("Residual weighting:", do_residual_weights), pos = 4, cex = CX)
            text(0, 0.025, paste0("Standardize RMSE:", do_standard), pos = 4, cex = CX)
          }
          next
        } else {
          if (set == "base") {
            para <- params_complete[1, ]
          } else {
            best <- get(set)
            para <- best$par
          }
          plot_diff <- refit(para, ml, square = sqorder)
          plot_title <- paste(set, "\n", paste0(nparams, ": ", round(denorma(para, nor), 2), collapse = " "))
        }

        # print(paste(set, round(mean(plot_diff), 2)))

        # to recompose the linear fit
        count <- count + 1

        im1 <- plot.prepare(lon, lat, plot_diff * fp$factor,
          proj = map_projection, lat_lim = lat_lim
        )
        im2 <- plot.prepare(lon, lat, plot_full * fp$factor,
          proj = map_projection, lat_lim = lat_lim
        )
        filled.contour3(im1$x, im1$y, im1$z,
          xlab = im1$xlab, ylab = im1$ylab, main = plot_title,
          levels = fp$lev_diff, color.palette = fp$color_diff, asp = im1$asp, extend = T,
          xlim = im1$xlim, ylim = im1$ylim, axes = im1$axes
        )
        contour(im2$x, im2$y, im2$z,
          levels = fp$lev_field,
          add = T, lwd = 2, drawlabels = T, labcex = clab
        )
        contour(im2$x, im2$y, im2$z,
          levels = -fp$lev_field,
          add = T, lwd = 2, drawlabels = T, lty = 2, labcex = clab
        )
        kol <- c("darkred", "darkgreen", "darkblue")
        for (region in regions) {
          rr <- sel_region(region)
          # bias <- weighted.mean(plot_diff, ww)
          rmse <- sqrt(weighted.sum(plot_diff[rr$lons, rr$lats]^2, area_ww[rr$lons, rr$lats]))
          text(120, 90 - which(region == regions) * 10, paste0(rr$sname, "RMSE:", round(rmse, 2)),
            cex = 3, col = kol[which(region == regions)]
          )
        }

        # text(120, 70, paste0("RMSE:", round(rmse, 2)), cex = 3, col = "darkred")
        # text(120, 80, paste0("NA RMSE:", round(rmse_na, 2)), cex = 3, col = "darkgreen")
        proj.addland(lon, lat, proj = map_projection)
        mtext(lettering[count], line = 1, adj = -0.05, cex = cex.letter)
        if (set == sets[length(sets) - 1]) {
          image.scale3(volcano,
            levels = fp$lev_diff, color.palette = fp$color_diff,
            colorbar.label = fp$legend_unit,
            cex.colorbar = imgscl_colorbar, cex.label = imgscl_label,
            colorbar.width = 1.5, line.colorbar = 0,
            line.label = fp$legend_distance
          )
        }
      }
    }

    dev.off()
  }
}

if (length(evolution$nens) > 8) {
  theplot <- list()
  evolution$grouper <- interaction(evolution$nens, evolution$region)
  for (nparam in nparams) {
    count <- which(nparam == nparams)
    theplot[[count]] <- ggplot(evolution, aes_string(x = "nens", y = nparam, fill = "region")) +
      geom_boxplot(aes(group = grouper)) +
      labs(
        title = paste(nparam, "evolution with increasing members"),
        x = paste("# members"), y = nparam, col = "Region"
      ) +
      stat_summary(mapping = aes(col = region), fun = mean, geom = "line", size = 2) +
      stat_summary(mapping = aes(label = sprintf("%0.2f", ..y..)), fun = mean, geom = "text", vjust = -0.25, size = 6) +
      coord_cartesian(xlim = c(10, 20), ylim = c(denorma(param_lower[nparam], nor) , 
                                                 denorma(param_upper[nparam], nor))) +
      theme_light(base_size = 18)

    if (count == length(nparams)) {
      plot_legend <- get_legend(theplot[[count]])
    }
    theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
  }
  name <- paste(FIGDIR, "/",
    "Spread_Evolution_params_ensemble", region, ".pdf",
    sep = ""
  )
  prow <- plot_grid(
    plotlist = theplot, align = "vh",
    labels = lettering[1:4], hjust = 0, nrow = 2,
    ncol = 2, label_size = 24
  )
  prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 2), nrow = 1, ncol = 2)
  save_plot(name, prow2, base_width = 20, base_height = 15)

  theplot <- list()
  procks <- c("rmse", "rsquared")
  for (prock in procks) {
    count <- which(prock == procks)
    if (prock == "rmse") {
      YY <- c(3, 7)
    } else {
      YY <- c(0, 0.5)
    }
    theplot[[count]] <- ggplot(evolution, aes_string(x = "nens", y = prock, fill = "region")) +
      geom_boxplot(aes(group = grouper)) +
      stat_summary(mapping = aes(label = sprintf("%0.3f", ..y..)), fun = mean, geom = "text", vjust = -0.25, size = 6) +
      labs(
        title = paste(prock, "evolution with increasing members"),
        x = paste("# members"), y = prock, col = "Region"
      ) +
      coord_cartesian(xlim = c(10, 20), ylim = YY) +
      theme_light(base_size = 18)

    if (count == length(procks)) {
      plot_legend <- get_legend(theplot[[count]])
    }
    theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
  }
  name <- paste(FIGDIR, "/",
    "Spread_Evolution_rmse_rsquared_ensemble", region, "_", season, ".pdf",
    sep = ""
  )
  prow <- plot_grid(
    plotlist = theplot, align = "vh",
    labels = lettering[1:4], hjust = 0, nrow = 2,
    ncol = 1, label_size = 24
  )
  prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 2), nrow = 1, ncol = 2)
  save_plot(name, prow2, base_width = 20, base_height = 12)
}

# print("ciao")
# theplot <- list()
# for (nparam in nparams) {
#  count <- which(nparam == nparams)
#  theplot[[count]] <- ggplot(evolution) +
#    geom_path(aes_string(x = "nens", y = nparam, col = "region"), size = 2) +
#    labs(
#      title = paste(nparam, "evolution with increasing members"),
#     x = paste("# members"), y = nparam, col = "Region"
#    ) +
#    coord_cartesian(xlim = c(3, 20), ylim = c(param_lower[nparam], param_upper[nparam])) +
#    theme_light(base_size = 18)

# if (count == length(nparams)) {
#   plot_legend <- get_legend(theplot[[count]])
# }
# theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
# }
# name <- paste(FIGDIR, "/", season, "/",
#  "Evolution_params_ensemble", region, "_", season, ".pdf",
#  sep = ""
# )
# prow <- plot_grid(
#  plotlist = theplot, align = "vh",
#  labels = lettering[1:4], hjust = 0, nrow = 2,
#  ncol = 2, label_size = 24
# )
# prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 2), nrow = 1, ncol = 2)
# save_plot(name, prow2, base_width = 20, base_height = 15)

# theplot <- list()
# procks <- c("rmse", "rsquared")
# for (prock in procks) {
#  count <- which(prock == procks)
#  if (prock == "rmse") {
#      YY <- c(0, 10)
#  } else {
#    YY <- c(0, 1)
#  }
#  theplot[[count]] <- ggplot(evolution) +
#    geom_path(aes_string(x = "nens", y = prock, col = "region"), size = 2) +
#    labs(
#      title = paste(prock, "evolution with increasing members"),
#      x = paste("# members"), y = prock, col = "Region"
#    ) +
#     coord_cartesian(xlim = c(3, 20), ylim = YY) +
#    theme_light(base_size = 18)
#
#    if (count == length(procks)) {
#      plot_legend <- get_legend(theplot[[count]])
#    }
#    theplot[[count]] <- theplot[[count]] + theme(legend.position = "none")
#  }
#  name <- paste(FIGDIR, "/", season, "/",
#    "Evolution_rmse_rsquared_ensemble", region, "_", season, ".pdf",
#  sep = ""
#  )
#  prow <- plot_grid(
#    plotlist = theplot, align = "vh",
#    labels = lettering[1:4], hjust = 0, nrow = 2,
#    ncol = 1, label_size = 24
#  )
#  prow2 <- plot_grid(prow, plot_legend, rel_widths = c(10, 2), nrow = 1, ncol = 2)
#  save_plot(name, prow2, base_width = 20, base_height = 12)
