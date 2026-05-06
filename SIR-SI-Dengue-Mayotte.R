# ============================================================================
# SIR-SI MODEL WITH CLIMATIC FORCING
# Comparison of 3 models : Constant, Lambrechts, Gaussian
# Available methods : MLE (Negative Binomial) ou WLS(Weighted Least Square)
# Applied method : MLE
# ============================================================================

# 1. Importing necessary packages

pkgs <- c("dplyr", "ggplot2", "tidyr", "deSolve", "gridExtra", 
          "readr", "lubridate", "scales", "knitr", "numDeriv")

new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new) > 0) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

# 2. Loading data
weekly_mayotte <- read.csv("data_complet_week.csv")
weekly_cases <- weekly_mayotte$n_cases
weekly_temp  <- weekly_mayotte$temp_mean_week
weekly_rain  <- weekly_mayotte$rain_mean_week
n_weeks <- length(weekly_cases)
times   <- seq(0, n_weeks * 7, by = 1)
l <- 7


# Fixed parameters
Nh <- 269579                    # Human population in 2019 (INSEE)
gamma_h <- 1/7                  # Human healing rate (1/7 per day)
mu_v <- 0.02                    # Mortality rate of vectors (1/45 per day)
k = 2                           # Average number of vector per human

# 3. Climatic forcing function
# 3.1 Lambrechts' Function (2011)
p_Lambrechts <- function(T) {
  out <- 0.001044 * T * (T - 12.286) * sqrt(pmax(32.461 - T, 0))
  out[out < 0 | !is.finite(out)] <- 0
  out
}

s_Lambrechts <- function(T) {
  raw <- (p_Lambrechts(T))^2
  max_raw <- max((p_Lambrechts(seq(15, 32, by = 0.1)))^2, na.rm = TRUE)
  out <- raw / max_raw
  out[!is.finite(out)] <- 0
  pmin(pmax(out, 0), 1)
}

# 3.2 Gaussian function with precipitation and temperature
g_temp <- function(T, T_opt, sigma_T) {
  exp(-(T - T_opt)^2 / (2 * sigma_T^2))
}

h_rain <- function(P, alpha_R) {
  1 + alpha_R * P
}

scaling_factor <- function(week, T_week, P_week, T_opt, sigma_T, alpha_R) {
  w_lag <- pmax(1, pmin(week - l, length(T_week)))
  Tt <- T_week[w_lag]
  Pt <- P_week[w_lag]
  out <- g_temp(Tt, T_opt, sigma_T) * h_rain(Pt, alpha_R)
  pmin(pmax(out, 0), 1)
}

# 4. SIR-SI model (5 COMPARTMENTS) =========================================
# dSh/dt = -beta(t) * Sh * Iv / Nh
# dIh/dt = beta(t) * Sh * Iv / Nh - gamma_h * Ih
# dRh/dt = gamma_h * Ih
# dSv/dt = mu_v * k * Nh - beta(t) * Sv * Ih / Nh - mu_v * Sv
# dIv/dt = beta(t) * Sv * Ih / Nh - mu_v * Iv

sir_si_ode <- function(t, y, pars) {
  Sh <- y[1]  # susceptible humans
  Ih <- y[2]  # Infected humans
  Rh <- y[3]  # Recovered humans
  Sv <- y[4]  # Susceptible vectors
  Iv <- y[5]  # Infected vectors
  
  # Calculation of beta(t) according to the model
  week <- pmax(1, pmin(ceiling(t / 7), pars$n_weeks))
  if(pars$model_type == "constant") {
    beta_t <- pars$beta0
  } else if(pars$model_type == "lambrechts") {
    beta_t <- pars$beta0 * s_Lambrechts(pars$temp_week[week])
  } else if(pars$model_type == "gaussian") {
    scale <- scaling_factor(week, pars$temp_week, pars$rain_week,
                            pars$T_opt, pars$sigma_T, pars$alpha_R)
    beta_t <- pars$beta0 * scale
  }
  
  beta_t <- max(beta_t, 1e-10)
  
  # Differential equations
  dSh <- -beta_t * Sh * Iv / (pars$Nh)
  dIh <- beta_t * Sh * Iv / (pars$Nh) - pars$gamma_h * Ih
  dRh <- pars$gamma_h * Ih
  dSv <- pars$mu_v * pars$k * pars$Nh - beta_t * Sv * Ih / pars$Nh - pars$mu_v * Sv
  dIv <- beta_t * Sv * Ih / pars$Nh - pars$mu_v * Iv
  
  list(c(dSh, dIh, dRh, dSv, dIv))
}

# 5. Predictions function
# calculation of the weekly incidence (new human cases)
incidence_weekly <- function(out) {
  daily_inc <- c(0, -diff(out$Sh))
  daily_inc[daily_inc < 0] <- 0
  w <- pmax(1, ceiling(out$time / 7))
  w <- pmin(w, length(weekly_cases))
  full <- tapply(daily_inc, w, sum, default = 0)
  full_out <- numeric(length(weekly_cases))
  full_out[seq_along(full)] <- full
  full_out
}

# Solve the model and return the predicted incidence
predict_incidence <- function(pars) {
  # Conditions initiales
  Iv0 <- pars$I0 / (pars$k * pars$Nh)  # Intitial infected vector
  y0 <- c(
    Sh0 = pars$Nh - pars$I0,
    Ih0 = pars$I0,
    Rh0 = 0,
    Sv0 = pars$k * pars$Nh - Iv0,
    Iv0 = Iv0
  )
  
  out <- tryCatch({
    ode(y = y0, times = times, func = sir_si_ode, 
        parms = pars, method = "lsoda", atol = 1e-6, rtol = 1e-6) |>
      as.data.frame()
  }, error = function(e) return(NULL))
  
  if(is.null(out)) return(rep(NA, length(weekly_cases)))
  
  names(out) <- c("time", "Sh", "Ih", "Rh", "Sv", "Iv")
  inc_week <- incidence_weekly(out)
  
  # Observed Incidence  (with observation factor rho)
  mu <- pmax(pars$rho * inc_week, 1e-9)
  mu[1:length(weekly_cases)]
}

# 6. Loss functions
# 6.1 MLE
neg_log_lik <- function(params_vec, model_type, T_opt = NULL, sigma_T = NULL) {
  # Parameters extraction
  beta0 <- params_vec[1]
  I0 <- params_vec[2]
  rho <- params_vec[3]
  phi <- params_vec[4]
  
  # Constraints
  if(beta0 <= 0.01 || beta0 > 2) return(1e10)
  if(I0 < 1 || I0 > 300) return(1e10)
  if(rho < 0.01 || rho > 0.3) return(1e10)
  if(phi <= 0.1 || phi > 100) return(1e10)
  
  # alpha_R parameter for the gaussian model
  alpha_R <- 0
  if(model_type == "gaussian") {
    if(length(params_vec) >= 5) {
      alpha_R <- params_vec[5]
      if(alpha_R < 0 || alpha_R > 10) return(1e10)
    }
    if(is.null(T_opt)) T_opt <- 29.0
    if(is.null(sigma_T)) sigma_T <- 4.0
  }
  
  # Construction of parameters model
  pars <- list(
    beta0 = beta0,
    gamma_h = gamma_h,
    mu_v = mu_v,
    k = k,
    Nh = Nh,
    I0 = I0,
    rho = rho,
    model_type = model_type,
    n_weeks = n_weeks,
    temp_week = weekly_temp,
    rain_week = weekly_rain
  )
  
  if(model_type == "gaussian") {
    pars$T_opt <- T_opt
    pars$sigma_T <- sigma_T
    pars$alpha_R <- alpha_R
  }
  
  # Prédiction
  mu <- predict_incidence(pars)
  if(any(is.na(mu))) return(1e10)
  
  # Log-likelihood Negative Binomial
  log_lik <- sum(dnbinom(weekly_cases, size = phi, mu = mu, log = TRUE))
  -log_lik
}


# 6.2 WLS
wls_loss <- function(params_vec, model_type, T_opt = NULL, sigma_T = NULL, 
                         initial_mu = NULL) {
  
  # Parameters extraction
  beta0 <- params_vec[1]
  I0 <- params_vec[2]
  rho <- params_vec[3]
  
  # Constraints
  if(beta0 <= 0.01 || beta0 > 2) return(1e10)
  if(I0 < 1 || I0 > 300) return(1e10)
  if(rho < 0.01 || rho > 0.3) return(1e10)
  
  alpha_R <- 0
  if(model_type == "gaussian") {
    if(length(params_vec) >= 4) {
      alpha_R <- params_vec[4]
      if(alpha_R < 0 || alpha_R > 10) return(1e10)
    }
    if(is.null(T_opt)) T_opt <- 29.0
    if(is.null(sigma_T)) sigma_T <- 4.0
  }
  
  # Model parameters
  pars <- list(
    beta0 = beta0,
    gamma_h = gamma_h,
    mu_v = mu_v,
    k = k,
    Nh = Nh,
    I0 = I0,
    rho = rho,
    model_type = model_type,
    n_weeks = n_weeks,
    temp_week = weekly_temp,
    rain_week = weekly_rain
  )
  
  if(model_type == "gaussian") {
    pars$T_opt <- T_opt
    pars$sigma_T <- sigma_T
    pars$alpha_R <- alpha_R
  }
  
  # Prediction
  mu <- predict_incidence(pars)
  if(any(is.na(mu))) return(1e10)
  
  # Observations
  obs <- weekly_cases
  
  # Weighted calculation
  weights <- 1 / pmax(mu, 0.5)
  
  # sum of weighted residual values
  wls <- sum(weights * (obs - mu)^2, na.rm = TRUE)
  
  return(wls)
}

# 6.3 Fonction générique qui choisit la méthode
choice_mle_wls <- function(params_vec, model_type, method = "mle", T_opt = NULL, sigma_T = NULL) {
  
  if(method == "mle") {
    return(neg_log_lik(params_vec, model_type, T_opt, sigma_T))
  } else if(method == "wls") {
    return(wls_loss(params_vec, model_type, T_opt, sigma_T))
  } else {
    stop("Choice 'mle' or 'wls'")
  }
}

# 7. Optimisation function
fit_model <- function(model_type, method = "mle", initial_params = NULL, 
                      T_opt = NULL, sigma_T = NULL, trace = 0) {
  
  cat("\n", rep("-", 50), "\n")
  cat("Model:", toupper(model_type), "- Method:", toupper(method), "\n")
  cat(rep("-", 50), "\n")
  
  # Initial parameters according to the model and the method
  if(is.null(initial_params)) {
    if(model_type == "gaussian") {
      if(method == "mle") {
        initial_params <- c(beta0 = 0.4, I0 = 10, rho = 0.1, phi = 1.0, alpha_R = 0.05)
      } else { # wls
        initial_params <- c(beta0 = 0.4, I0 = 10, rho = 0.1, alpha_R = 0.05)
      }
    } else {
      if(method == "mle") {
        initial_params <- c(beta0 = 0.4, I0 = 10, rho = 0.1, phi = 1.0)
      } else { # wls
        initial_params <- c(beta0 = 0.4, I0 = 10, rho = 0.1)
      }
    }
  }
  
  # Bounds according to the model and the method
  if(model_type == "gaussian") {
    if(method == "mle") {
      lower <- c(beta0 = 0.01, I0 = 1, rho = 0.01, phi = 0.1, alpha_R = 0)
      upper <- c(beta0 = 2, I0 = 300, rho = 0.3, phi = 30, alpha_R = 10)
    } else { # wls
      lower <- c(beta0 = 0.01, I0 = 1, rho = 0.01, alpha_R = 0)
      upper <- c(beta0 = 2, I0 = 300, rho = 0.3, alpha_R = 10)
    }
  } else {
    if(method == "mle") {
      lower <- c(beta0 = 0.01, I0 = 1, rho = 0.01, phi = 0.1)
      upper <- c(beta0 = 2, I0 = 300, rho = 0.3, phi = 30)
    } else { # wls
      lower <- c(beta0 = 0.01, I0 = 1, rho = 0.01)
      upper <- c(beta0 = 2, I0 = 300, rho = 0.3)
    }
  }
  
  # Optimisation
  fit <- optim(
    par = initial_params,
    fn = function(p) choice_mle_wls(p, model_type = model_type, method = method, 
                                   T_opt = T_opt, sigma_T = sigma_T),
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 500, trace = trace)
  )
  
  if(method == "mle") {
    fit$logLik <- -fit$value
    fit$AIC <- 2 * length(fit$par) - 2 * fit$logLik
    fit$BIC <- log(n_weeks) * length(fit$par) - 2 * fit$logLik
    fit$metric_name <- "Log-Likelihood"
    fit$metric_value <- fit$logLik
  } else {
    fit$WLS <- fit$value
    fit$RSE <- sqrt(fit$WLS / (n_weeks - length(fit$par)))
    fit$metric_name <- "wls"
    fit$metric_value <- fit$WLS
  }
  
  fit$model_type <- model_type
  fit$method <- method
  
  return(fit)
}

# 8.  Implementation of optimisations
# Method (= mle ou wls)
# The choice is mle here
fit_constant <- fit_model("constant", method = "mle")
fit_lambrechts <- fit_model("lambrechts", method = "mle")
fit_gaussian <- fit_model("gaussian", method = "mle", T_opt = 29.0, sigma_T = 4.0)


# 9. Prediction
make_prediction <- function(fit) {
  pars <- list(
    model_type = fit$model_type,
    beta0 = fit$par["beta0"],
    I0 = fit$par["I0"],
    rho = fit$par["rho"],
    gamma_h = gamma_h,
    mu_v = mu_v,
    k = k, 
    Nh = Nh,
    n_weeks = n_weeks,
    temp_week = weekly_temp,
    rain_week = weekly_rain
  )
  
  if(fit$model_type == "gaussian") {
    pars$T_opt <- 29.0
    pars$sigma_T <- 4.0
    pars$alpha_R <- fit$par["alpha_R"]
  }
  
  predict_incidence(pars)
}

# Constant model
pred_constant <- make_prediction(fit_constant)
# Lambrechts' model
pred_lambrechts <- make_prediction(fit_lambrechts)
# Gaussian model with climatic forcing (temperature and precipitation)
pred_gaussian <- make_prediction(fit_gaussian)

# 10. Tables of results
params <- data.frame(
  Model = c("Constant", "Lambrechts", "Gaussien"),
  beta0 = c(round(fit_constant$par["beta0"], 3),
            round(fit_lambrechts$par["beta0"], 3),
            round(fit_gaussian$par["beta0"], 3)),
  I0 = c(round(fit_constant$par["I0"], 1),
         round(fit_lambrechts$par["I0"], 1),
         round(fit_gaussian$par["I0"], 1)),
  rho = c(round(fit_constant$par["rho"], 3),
          round(fit_lambrechts$par["rho"], 3),
          round(fit_gaussian$par["rho"], 3)),
  phi = c(round(fit_constant$par["phi"], 2),
          round(fit_lambrechts$par["phi"], 2),
          round(fit_gaussian$par["phi"], 2)),
  alpha_R = c("-", "-", round(fit_gaussian$par["alpha_R"], 3)),
  LogLik = c(round(fit_constant$logLik, 1),
             round(fit_lambrechts$logLik, 1),
             round(fit_gaussian$logLik, 1)),
  AIC = c(round(fit_constant$AIC, 1),
          round(fit_lambrechts$AIC, 1),
          round(fit_gaussian$AIC, 1))
)

# 11. Visualisations
df_plot <- data.frame(
  Week = 1:n_weeks,
  Observed = weekly_cases,
  Constant = pred_constant,
  Lambrechts = pred_lambrechts,
  Gaussian = pred_gaussian
)

df_long <- pivot_longer(df_plot, cols = c("Constant", "Lambrechts", "Gaussian"),
                            names_to = "Model", values_to = "Predict")

p <- ggplot() +
  geom_col(data = df_plot, aes(x = Week, y = Observed), 
           fill = "#FFBEB2", alpha = 0.6, width = 0.8) +
  geom_line(data = df_long, aes(x = Week, y = Predict, color = Model), size = 1.2) +
  scale_color_manual(values = c("Constant" = "#009E73", "Lambrechts" = "#AE123A", "Gaussian" = "#0072B2")) +
  labs(x = "Week", y = "Incidence", title = "SIR-SI")+ 
  theme(
      legend.position="bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )

print(p)

cat("Model selected by AIC: ", params$Model[which.min(params$AIC)], "\n")


# 12. Visualisations of lag and basic reproduction number (Through beta(t))
# normalize variables between 0 and 1
weekly_norm <- weekly_mayotte %>%
  mutate(
    cases_norm = scales::rescale(n_cases),
    temp_norm  = scales::rescale(temp_mean_week),
    hum_norm   = scales::rescale(hum_mean_week),
    rain_norm  = scales::rescale(rain_mean_week)
  )

ts_cases <- weekly_norm$cases_norm
ts_temp  <- weekly_norm$temp_norm
ts_hum   <- weekly_norm$hum_norm
ts_rain  <- weekly_norm$rain_norm

# compute correlations for lags 0 to 12
lags <- 0:30
lag_corr <- tibble(
  lag = lags,
  temp_corr = sapply(lags, function(l) cor(ts_cases[(l+1):length(ts_cases)],
                                           ts_temp[1:(length(ts_temp)-l)],
                                           use="complete.obs")),
  hum_corr  = sapply(lags, function(l) cor(ts_cases[(l+1):length(ts_cases)],
                                           ts_hum[1:(length(ts_hum)-l)],
                                           use="complete.obs")),
  rain_corr = sapply(lags, function(l) cor(ts_cases[(l+1):length(ts_cases)],
                                           ts_rain[1:(length(ts_rain)-l)],
                                           use="complete.obs"))
) %>%
  pivot_longer(cols = -lag, names_to = "variable", values_to = "correlation")

print(ggplot(lag_corr, aes(x = lag, y = correlation, color = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "temp_corr" = "#AE123A",
    "hum_corr" = "#FFC685",
    "rain_corr" = "#0072B2"
  )) +
  labs(
    title = "Correlation Between Climate Variables and Dengue Cases by Lag (0–12 weeks)",
    x = "Lag (weeks)",
    y = "Correlation coefficient",
    color = "Variable"
  )  +
  theme(
    legend.position="bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ))


# Parameters
beta0 <- params$beta0[3]
alpha_R <- params$beta0[3]
l <- 7
T_opt <- 29
sigma_T <- 4
gamma <- 1/7
mu_v <- 0.02
k <- 2

# Threshold
seuil <- sqrt((gamma * mu_v)/k)

weeks <- 1:n_weeks
T_lag <- c(rep(NA, l), weekly_temp[1:(n_weeks - l)])
P_lag <- c(rep(NA, l), weekly_rain[1:(n_weeks - l)])

# Climatic functions
g <- function(T) {
  exp(-((T - T_opt)^2) / (2 * sigma_T^2))
}

h <- function(P) {
  1 + alpha_R * P
}

# beta(t)
beta_t <- beta0 * g(T_lag) * h(P_lag)
df <- data.frame(
  week = weeks,
  beta_t = beta_t,
  threshold = seuil
) %>%
  na.omit() %>%
  mutate(
    regime = ifelse(beta_t > threshold, "Above threshold", "Below threshold")
  )

# Graphics
p <- ggplot(df, aes(x = week, y = beta_t)) +
  geom_ribbon(
    data = subset(df, beta_t > threshold),
    aes(ymin = threshold, ymax = beta_t),
    fill = "red", alpha = 0.2
  ) +
  geom_line(aes(color = regime), linewidth = 1) +
  geom_hline(
    yintercept = seuil,
    linetype = "dashed",
    linewidth = 1.2,
    color = "black"
  ) +
  scale_color_manual(values = c(
    "Above threshold" = "red",
    "Below threshold" = "blue"
  )) +
  labs(
    title = expression("Dynamics of " * beta(t)),
    x = "Week (March 2019 - July 2020)",
    y = expression(beta(t)),
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)
