# source('script/data_preprocessing_new.R')
rm(list = ls())
suppressMessages({source('PIV.R')}); set.seed(2026); 
load("results/ALS_DATA.RData"); 
truncate <- function(x){pmax(pmin(x, quantile(x, 0.95)), quantile(x, 0.05))}
DATA <- DATA %>% mutate(
  year_edu  = as.numeric(year_edu > 12),
  wage_hour = pmin(wage_hour, 30),
  age       = truncate(age),
  # age       = as.numeric(age > 20),
  attitude  = truncate(attitude),
  year_expe = truncate(year_expe),
); par(mfrow=c(2,3))
hist(DATA$wage_hour,breaks = 10,main='Hourly wages')
hist(DATA$age, main= 'Age', breaks = 10)
hist(DATA$attitude, main= 'Attitude', breaks = 10)
hist(log(DATA$year_expe), main= 'Experience', breaks = 10)
hist(DATA$year_edu, main= 'year_edu', breaks = 10)
hist(DATA$year_expe, main= 'year_expe', breaks = 10)
head(DATA);summary(DATA);summary(lm(wage_hour ~ year_edu, data = DATA))


set.seed(123)
formula = "wage_hour ~ year_edu ~ sex | t + attitude + capital + born_australia"
# formula = "wage_hour ~ year_edu ~ sex | attitude + capital + born_australia"
# formula = "wage_hour ~ year_edu ~ sex | t + attitude + born_australia"
# formula = "wage_hour ~ year_edu ~ age + sex | t + attitude"
fit.PIV <- PIV_mgcv(
  formula = formula,data = DATA, invalid = c(0,1,2), 
  print_state = 0, K = 1, smoothness = 5
); print(fit.PIV) # very good
fit.GE0 <- G_estimation(formula = formula, data = DATA, invalid = 2, smoothness = 5);print(fit.GE0)
fit.GE1 <- G_estimation(formula = formula, data = DATA, invalid = 1, smoothness = 5);print(fit.GE1)
fit.GE2 <- G_estimation(formula = formula, data = DATA, invalid = 0, smoothness = 5);print(fit.GE2)

set.seed(123)
DML <- DoubleML::DoubleMLIRM$new(
  data = DoubleML::double_ml_data_from_data_frame(
    df = as.data.frame(DATA),y_col = "wage_hour",d_cols = "year_edu",
    x_cols = c("t", "attitude", "capital", "born_australia",'sex')
  ),
  ml_g = mlr3::lrn("regr.ranger",num.trees = 100,min.node.size = 2,max.depth = 5),
  ml_m = mlr3::lrn("classif.ranger",num.trees = 100,min.node.size = 2,max.depth = 5)
)
DML$fit()
fit.DML <- DML$summary()







results <- data.frame(
  PIV4 = c(est = fit.PIV$res$est_aug[1], sd = fit.PIV$res$sd_aug[1]),
  PIV3 = c(est = fit.PIV$res$est_aug[2], sd = fit.PIV$res$sd_aug[2]),
  PIV2 = c(est = fit.PIV$res$est_aug[3], sd = fit.PIV$res$sd_aug[3]),
  GE4 = c(est = fit.GE2[[1]], sd = fit.GE2[[2]]),
  GE3 = c(est = fit.GE1[[1]], sd = fit.GE1[[2]]),
  GE2 = c(est = fit.GE0[[1]], sd = fit.GE0[[2]]),
  NUC = c(est = fit.DML[1], sd = fit.DML[2])
) %>%
  t() %>%
  as.data.frame()

method_order <- c("PIV2", "PIV3", "PIV4", "GE2", "GE3", "GE4", "NUC")
results <- results[method_order, , drop = FALSE] %>%
  mutate(
    p_value = 2 * pnorm(-abs(est / sd)),
    ci_lower = est - qnorm(0.975) * sd,
    ci_upper = est + qnorm(0.975) * sd
  )
print(results)


############################
## LaTeX results table
############################
format_fixed <- function(x, digits) {
  formatC(x, format = "f", digits = digits)
}

latex_table <- paste0(
  "\\begin{table}[t]\n",
  "    \\centering\n",
  "    \\begin{tabular}{c|ccc ccc c}\n",
  "        \\toprule\n",
  "        Method & PIV$_2$ & PIV$_3$ & PIV$_4$ & GE$_2$ & GE$_3$ & GE$_4$ & NUC \\\\\n",
  "        \\midrule\n",
  "        EST & ", paste(format_fixed(results$est, 2), collapse = " & "), " \\\\\n",
  "        SE & ", paste(format_fixed(results$sd, 2), collapse = " & "), " \\\\\n",
  "        $p$-value & ", paste(ifelse(
    format_fixed(results$p_value, 3) == '0.000',
    '<0.001',
    format_fixed(results$p_value, 3)
  ), collapse = " & "), " \\\\\n",
  "        \\bottomrule\n",
  "    \\end{tabular}\n",
  "    \\caption{Estimated effects of education on hourly wages.}\n",
  "    \\label{tbl: ALS_results}\n",
  "\\end{table}\n"
)

cat(latex_table)
writeLines(latex_table, "results/ALS_results_table.tex")


############################################
## Point estimates and confidence intervals
############################################
plot_results <- results %>%
  tibble::rownames_to_column("Method") %>%
  mutate(Method = factor(Method, levels = method_order))

estimate_plot <- ggplot2::ggplot(
  plot_results,
  ggplot2::aes(x = Method, y = est, ymin = ci_lower, ymax = ci_upper)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbar(width = 0.2, linewidth = 0.3) +
  ggplot2::geom_point(size = 1) +
  ggplot2::labs(
    x = NULL,
    y = "Estimated effect on hourly wage"
  ) +
  ggplot2::theme_classic(base_size = 12)

print(estimate_plot)
ggplot2::ggsave(
  filename = "results/ALS_estimates_CI.pdf",
  plot = estimate_plot,
  width = 7,
  height = 3
)



