
library(dplyr)
library(ggplot2)
library(purrr)

for(sample in c(2000,5000,10000)){
  ALL_RES <- list()
  for(i in 1:6){
    
    # path <- paste0("data/3_10000/RES", i, ".rds")
    path <- paste0("results/", sample, "/RES", i, "_5.rds")
    RES <- readRDS(path)

    RES$DGP <- paste0("DGP ", i)

    ALL_RES[[i]] <- RES
  }

  plot_dat <- bind_rows(ALL_RES)

  plot_dat$method <- factor(
    plot_dat$method,
    levels = c(
      "PIV1","PIV2","PIV3",
      "GE1","GE2","GE3",
      "NUC","ivreg","TSHS","sisVIVE"
    ),
    labels = c(
      "PIV1","PIV2","PIV3",
      "GE1","GE2","GE3",
      "NUC","2SLS","TSHS","sisVIVE"
    )
  )

  ## =========================
  ## 5% trimming
  ## =========================
  plot_dat_trim <- plot_dat %>%
    group_by(DGP, method) %>%
    filter(
      est >= quantile(est, 0.05, na.rm = TRUE),
      est <= quantile(est, 0.95, na.rm = TRUE)
    ) %>%
    ungroup()

  ## =========================
  ## violin plot
  ## =========================
  p <- ggplot(
    plot_dat_trim,
    aes(x = method, y = est)
  ) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(
      width = 0.12,
      outlier.size = 0.3,
      fill = "white"
    ) +
    facet_wrap(~DGP, nrow = 2, ncol = 3, scales = "free_y") +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      color = "red"
    ) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Violin plot of estimate bias (5% trimmed)"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = element_text(size = 12)
    )

  ggsave(
    filename = paste0("results/fig/violinplot_3_", sample, ".pdf"),
    plot = p,
    width = 12,
    height = 6
  )
}
