# big R function to loop on what we neeed
# clean everything to be sure
rm(list = ls())

# directory
TUNEDIR <- "/home/paolo/RScript/reforge/tune_gwd/aggressive_tuning"

# loop on var and levels
big_vars <- c("ua")
big_levels <- list(
  c(5000, 25000, 85000)
)

# loop on parameters
big_nparams <- list( # c("GKWAKE", "GKDRAG", "GFRCRIT"),
  c("ZTOFD", "GKWAKE", "GKDRAG", "GFRCRIT")
)

# loop on seasons
big_seasons <- list(
  "DJFM",
  c("DJFM", "JJAS")
)

# residual and standardization
big_residual <- c(FALSE)
big_standard <- c(FALSE, TRUE)

# linear or quadratic fit
big_order <- c("linear")

# parameters to be excluded
big_params_exclude <- list(
  NULL,
  "GFRCRIT",
  "GKDRAG"
)

#### NO LOOP ###

# regions
regions <- c("global", "northern_hemisphere")


# greater loops
for (vars in big_vars) {
  for (ff in seq_along(big_nparams)) {
    for (kk in seq_along(big_seasons)) {
      for (zz in seq_along(big_levels)) {
        for (xx in seq_along(big_params_exclude)) {
          levels <- big_levels[[zz]]
          seasons <- big_seasons[[kk]]
          nparams <- big_nparams[[ff]]
          params_exclude <- big_params_exclude[[xx]]
          for (do_residual_weights in big_residual) {
            for (do_standard in big_standard) {
              for (sqorder in big_order) {
                print(paste(
                  vars, paste(levels, collapse = "-"), paste(seasons, collapse = "-"),
                  paste(nparams, collapse = "-"), do_residual_weights, do_standard, sqorder
                ))
                source(file.path(TUNEDIR, "new_fingerprint.R"))
              }
            }
          }
        }
      }
    }
  }
}
