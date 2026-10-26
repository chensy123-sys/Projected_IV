rm(list = ls())

parse_integer_argument <- function(value, name){
  parsed <- suppressWarnings(as.numeric(value))
  if(length(parsed) != 1L || is.na(parsed) || !is.finite(parsed) ||
     parsed != floor(parsed) || abs(parsed) > .Machine$integer.max){
    stop(sprintf("%s must be an integer", name), call. = FALSE)
  }
  as.integer(parsed)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2L){
  stop(
    "Usage: Rscript script/monteCarlo.R <case> <n>",
    call. = FALSE
  )
}

case <- parse_integer_argument(args[[1]], "case")
n <- parse_integer_argument(args[[2]], "n")

if(!case %in% 1:6){
  stop("case must be one of 1, 2, 3, 4, 5, or 6", call. = FALSE)
}
if(n < 1L){
  stop("n must be a positive integer", call. = FALSE)
}

# Resolve paths relative to this script so it can be launched from any directory.
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if(length(file_arg) != 1L){
  stop("Unable to determine the path of monteCarlo.R", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", file_arg[[1]]),
  mustWork = TRUE
)
script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressMessages(source(file.path(project_dir, "PIV.R")))
source(file.path(script_dir, "info.R"))

results_dir <- file.path(project_dir, "results", sprintf("%d_%d", n, P))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
max_try <- 5

message(sprintf("Running DGP %d with n = %d", case, n))

for(loop in seq_len(Looptime)){
  for(attempt in seq_len(max_try)){
    ok <- TRUE

    tryCatch({
      sim_dat <- sim_PIV(n, case = case, P = P, ATE = ATE)

      invisible(suppressMessages(suppressWarnings({
        time2 <- system.time({
          # PIV
          fit.PIV <- PIV_mgcv(
            formula = "Y ~ A ~ X | Z1 + Z2 + Z3",
            data = sim_dat,
            invalid = c(0, 1, 2)
          )

          # GE
          fit.GE3 <- G_estimation(
            formula = Y ~ A ~ X | Z1 + Z2 + Z3,
            data = sim_dat,
            invalid = 0
          )
          fit.GE2 <- G_estimation(
            formula = Y ~ A ~ X | Z1 + Z2 + Z3,
            data = sim_dat,
            invalid = 1
          )
          fit.GE1 <- G_estimation(
            formula = Y ~ A ~ X | Z1 + Z2 + Z3,
            data = sim_dat,
            invalid = 2
          )

          # 2SLS
          if(case %in% c(1, 3, 5)){
            fit.ivreg <- summary(AER::ivreg(
              formula = Y ~ A + X + Z1 + Z2 | Z1 + Z2 + Z3 + X,
              data = sim_dat
            ))
          }else{
            fit.ivreg <- summary(AER::ivreg(
              formula = Y ~ A + X + Z1 | Z1 + Z2 + Z3 + X,
              data = sim_dat
            ))
          }

          # DML (NUC)
          DML <- DoubleML::DoubleMLIRM$new(
            data = DoubleML::double_ml_data_from_data_frame(
              df = as.data.frame(sim_dat),
              y_col = "Y",
              d_cols = "A",
              x_cols = c("Z1", "Z2", "Z3", "X")
            ),
            ml_g = mlr3::lrn(
              "regr.ranger",
              num.trees = 100,
              min.node.size = 2,
              max.depth = 5
            ),
            ml_m = mlr3::lrn(
              "classif.ranger",
              num.trees = 100,
              min.node.size = 2,
              max.depth = 5
            )
          )
          DML$fit()
          fit.DML <- DML$summary()

          # sisVIVE
          res_Y <- residuals(lm(Y ~ X, data = sim_dat))
          res_D <- residuals(lm(A ~ X, data = sim_dat))
          res_Z <- residuals(
            lm(cbind(Z1, Z2, Z3) ~ X, data = sim_dat)
          ) |> as.matrix()
          fit.sisVIVE <- sisVIVE::cv.sisVIVE(
            Y = res_Y,
            D = res_D,
            Z = res_Z,
            K = 3
          )

          # TSHT
          fit.TSHT <- summary(RobustIV::TSHT(
            Y = sim_dat$Y,
            D = sim_dat$A,
            Z = as.matrix(sim_dat %>% dplyr::select(Z1, Z2, Z3)),
            X = sim_dat$X
          ))
        })
      })))
    }, error = function(e){
      ok <<- FALSE
      message(sprintf(
        "loop=%d attempt=%d failed: %s",
        loop,
        attempt,
        e$message
      ))
      Sys.sleep(2)
    })

    if(ok){
      res_list <- list(
        PIV3 = c(est = fit.PIV$res$est_aug[1], sd = fit.PIV$res$sd_aug[1]),
        PIV2 = c(est = fit.PIV$res$est_aug[2], sd = fit.PIV$res$sd_aug[2]),
        PIV1 = c(est = fit.PIV$res$est_aug[3], sd = fit.PIV$res$sd_aug[3]),
        GE3 = c(est = fit.GE3[[1]], sd = fit.GE3[[2]]),
        GE2 = c(est = fit.GE2[[1]], sd = fit.GE2[[2]]),
        GE1 = c(est = fit.GE1[[1]], sd = fit.GE1[[2]]),
        NUC = c(est = fit.DML[1], sd = fit.DML[2]),
        ivreg = c(
          est = fit.ivreg$coefficients["A", "Estimate"],
          sd = fit.ivreg$coefficients["A", "Std. Error"]
        ),
        sisVIVE = c(est = fit.sisVIVE$beta, sd = NA),
        TSHS = c(est = fit.TSHT[1, 1], sd = fit.TSHT[1, 2])
      )

      res_all <- do.call(rbind, res_list) |> as.data.frame()
      res_all$method <- rownames(res_all)
      rownames(res_all) <- NULL
      res_all <- res_all[, c("method", "est", "sd")]
      res_all$est <- res_all$est - ATE
      res_df <- res_all |> dplyr::mutate(loop = loop)

      if(!exists("res")){
        res <- res_df
      }else{
        res <- rbind(res, res_df)
      }

      print(loop)
      break
    }

    if(attempt == max_try){
      message(sprintf(
        "loop=%d failed after %d attempts",
        loop,
        max_try
      ))
    }
  }
}

if(!exists("res")){
  stop("All Monte Carlo attempts failed; no result file was written")
}

RES <- data.frame(res, row.names = seq_len(nrow(res)))
output_file <- file.path(results_dir, sprintf("RES%d.rds", case))
saveRDS(RES, output_file)
message(sprintf("Saved results to %s", output_file))
