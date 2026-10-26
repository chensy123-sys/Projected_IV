rm(list = ls())

library(ggplot2)

set.seed(123)
n <- 500
U <- rnorm(n)
Z <- runif(n, min = -3, max = 3)
A <- Z^3 + U - 0.5 * Z
Y <- A + Z + U + 1

beta_values <- c(0.5, 0.75, 1, 1.25, 1.5)

plot_data <- do.call(
  rbind,
  lapply(beta_values, function(beta_value) {
    data.frame(
      Z = Z,
      residual = Y - beta_value * A,
      beta = factor(beta_value, levels = beta_values)
    )
  })
)

test_Z <- seq(-3, 3, by = 0.1)
line_data <- do.call(
  rbind,
  lapply(beta_values, function(beta_value) {
    fit_data <- data.frame(Z = Z, residual = Y - beta_value * A)
    fit <- lm(residual ~ Z, data = fit_data)

    data.frame(
      Z = test_Z,
      fitted = predict(fit, newdata = data.frame(Z = test_Z)),
      beta = factor(beta_value, levels = beta_values)
    )
  })
)

title_data <- data.frame(
  Z = 0,
  beta = factor(beta_values, levels = beta_values),
  label = paste0("beta[a] == ", beta_values),
  title_color = ifelse(beta_values == 1, "red", "#222222")
)

p <- ggplot(plot_data, aes(x = Z, y = residual)) +
  geom_point(color = "#4D4D4D", alpha = 0.5, size = 0.5) +
  geom_line(
    data = line_data,
    aes(y = fitted),
    color = "blue",
    linewidth = 1
  ) +
  geom_text(
    data = title_data,
    aes(x = Z, y = Inf, label = label, color = title_color),
    vjust = 1.4,
    fontface = "bold",
    size = 4.2,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~beta, nrow = 1, scales = "free_y") +
  scale_color_identity() +
  labs(
    x = "Z",
    y = expression(Y - beta[a] * A)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.border = element_rect(color = "#444444", fill = NA, linewidth = 0.6),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

print(p)

ggsave(
  filename = "results/linear_model.pdf",
  plot = p,
  width = 12,
  height = 2.5,
  dpi = 300
)
























rm(list = ls())

library(ggplot2)

set.seed(123)
n <- 500
U <- rnorm(n)
Z <- runif(n, min = -3, max = 3)
A <- Z^3 + U - 0.5 * Z
Y <- A + Z^2 + U + 1

beta_values <- c(0.5, 0.75, 1, 1.25, 1.5)

plot_data <- do.call(
  rbind,
  lapply(beta_values, function(beta_value) {
    data.frame(
      Z = Z,
      residual = Y - beta_value * A,
      beta = factor(beta_value, levels = beta_values)
    )
  })
)

test_Z <- seq(-3, 3, by = 0.1)
line_data <- do.call(
  rbind,
  lapply(beta_values, function(beta_value) {
    fit_data <- data.frame(Z = Z, residual = Y - beta_value * A)
    fit <- lm(residual ~ Z + I(Z^2), data = fit_data)
    
    data.frame(
      Z = test_Z,
      fitted = predict(fit, newdata = data.frame(Z = test_Z)),
      beta = factor(beta_value, levels = beta_values)
    )
  })
)

title_data <- data.frame(
  Z = 0,
  beta = factor(beta_values, levels = beta_values),
  label = paste0("beta[a] == ", beta_values),
  title_color = ifelse(beta_values == 1, "red", "#222222")
)

p <- ggplot(plot_data, aes(x = Z, y = residual)) +
  geom_point(color = "#4D4D4D", alpha = 0.5, size = 0.5) +
  geom_line(
    data = line_data,
    aes(y = fitted),
    color = "blue",
    linewidth = 1
  ) +
  geom_text(
    data = title_data,
    aes(x = Z, y = Inf, label = label, color = title_color),
    vjust = 1.4,
    fontface = "bold",
    size = 4.2,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~beta, nrow = 1, scales = "free_y") +
  scale_color_identity() +
  labs(
    x = "Z",
    y = expression(Y - beta[a] * A)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.border = element_rect(color = "#444444", fill = NA, linewidth = 0.6),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(8, 10, 8, 8)
  )

print(p)

ggsave(
  filename = "results/quadratic_model.pdf",
  plot = p,
  width = 12,
  height = 2.5,
  dpi = 300
)
