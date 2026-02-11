 # ============================================================
 # SeaTac temperature comparison: 1950 vs 2021
 # - Reads NOAA-style daily data from seatac_data.csv
 # - Plots daily MAX temp + sinusoidal fit for 1950 and 2021
 # - Plots daily MIN temp + sinusoidal fit for 1950 and 2021
 # - Plots daily MEAN temp (max+min)/2 + sinusoidal fit for 1950 and 2021
 #
 # Notes:
 # - NOAA daily temps in this file are in tenths of °C (common GHCN format).
 # - Convert to °F: (C * 9/5) + 32, where C = value/10
 # - Sinusoid fit uses: T(d) = a*sin(2πd/365) + b*cos(2πd/365) + c
 #   then converts to amplitude/phase form.
 # ============================================================

 library(tidyverse)
 library(lubridate)

 # ----------------------------
 # 1) Read in the data
 # ----------------------------
 file_path <- "seatac_data.csv"
#
 raw <- read_csv(file_path, show_col_types = FALSE)

 # Helper: tenths of °C -> °F
 to_F <- function(x_tenthsC) {
   (x_tenthsC / 10) * 9/5 + 32
 }

 # Parse dates + compute daily max/min/mean in °F
 # (Require both TMAX and TMIN for mean; allow TMAX-only for max plot, etc.)
 dat <- raw %>%
   mutate(
     DATE = mdy(DATE),
     year = year(DATE),
     doy  = yday(DATE),
     tmax_F = if_else(!is.na(TMAX), to_F(TMAX), NA_real_),
     tmin_F = if_else(!is.na(TMIN), to_F(TMIN), NA_real_),
     tmean_F = if_else(!is.na(TMAX) & !is.na(TMIN), to_F((TMAX + TMIN) / 2), NA_real_)
   ) %>%
   filter(year %in% c(1950, 2021))

 # ----------------------------
 # 2) Sinusoidal fit function
 # ----------------------------
 # Fit: T(d) = a*sin(w d) + b*cos(w d) + c, w = 2π/365
 # Return coefficients + a tidy prediction dataframe for plotting.
 fit_sinusoid <- function(df_year, value_col) {
   w <- 2 * pi / 365

   # Keep complete cases for the chosen variable
   d <- df_year %>%
     select(doy, value = all_of(value_col)) %>%
     filter(!is.na(value))

   # Linear regression in sin/cos basis
   fit <- lm(value ~ sin(w * doy) + cos(w * doy), data = d)

   coefs <- coef(fit)
   c0 <- unname(coefs[1])
   a  <- unname(coefs[2])
   b  <- unname(coefs[3])

   # Convert to amplitude/phase: a*sin + b*cos = A*sin(w d + phi)
 A   <- sqrt(a^2 + b^2)
   phi <- atan2(b, a)  # radians

   # Prediction grid across a full year
   grid <- tibble(doy = 1:365) %>%
     mutate(
       pred = c0 + a * sin(w * doy) + b * cos(w * doy)
     )

   list(
     fit = fit,
     params = tibble(A = A, phi = phi, C = c0),
     grid = grid
   )
 }

 # ----------------------------
 # 3) Plot helper (daily dots + fitted curve)
 # ----------------------------
 plot_daily_with_fit <- function(df, value_col, y_label, title_text) {
   # Fit per year
   fits <- df %>%
     group_split(year) %>%
     set_names(map_chr(., ~ as.character(unique(.x$year)))) %>%
     map(~ fit_sinusoid(.x, value_col))

   # Build prediction dataframe with year column
   pred_df <- imap_dfr(fits, ~ .x$grid %>% mutate(year = as.integer(.y)))

   # Optional: print equations to console (nice for copying into notes)
   # T(d) = A sin(2π d/365 + phi) + C
   eqs <- imap_dfr(fits, ~ {
     p <- .x$params
     tibble(
       year = as.integer(.y),
       A = p$A,
       phi = p$phi,
       C = p$C,
       equation = sprintf("T(d) = %.2f * sin(2π d/365 + %.2f) + %.2f", A, phi, C)
     )
   })
   print(eqs)

   ggplot(df, aes(x = doy, y = .data[[value_col]], group = factor(year))) +
     geom_point(aes(color = factor(year)), alpha = 0.35, size = 0.9) +
     geom_line(
       data = pred_df,
       aes(x = doy, y = pred, color = factor(year)),
       linewidth = 1.1
     ) +
     labs(
       x = "Day of Year",
       y = y_label,
       color = "Year",
       title = title_text
     ) +
     theme_minimal(base_size = 12)
 }

 # ----------------------------
 # 4) Plot DAILY MAX for 1950 and 2021 + fits
 # ----------------------------
 p_max <- plot_daily_with_fit(
   df = dat,
   value_col = "tmax_F",
   y_label = "Daily Maximum Temperature (°F)",
   title_text = "SeaTac Daily Maximum Temperature: 1950 vs 2021 (Sinusoidal Fits)"
 )
 print(p_max)

 # ----------------------------
 # 5) Plot DAILY MIN for 1950 and 2021 + fits
# # ----------------------------
 p_min <- plot_daily_with_fit(
   df = dat,
   value_col = "tmin_F",
   y_label = "Daily Minimum Temperature (°F)",
   title_text = "SeaTac Daily Minimum Temperature: 1950 vs 2021 (Sinusoidal Fits)"
 )
 print(p_min)

 # ----------------------------
 # 6) Plot DAILY MEAN ( (max+min)/2 ) for 1950 and 2021 + fits
 # ----------------------------
 p_mean <- plot_daily_with_fit(
   df = dat,
   value_col = "tmean_F",
   y_label = "Daily Mean Temperature (°F)",
   title_text = "SeaTac Daily Mean Temperature: 1950 vs 2021 (Sinusoidal Fits)"
 ) +
   geom_hline(
     yintercept = 65,
     linetype = "dashed",
     linewidth = 0.9,
     color = "black"
   )

 print(p_mean)

# ----------------------------
# (Optional) Save plots
# ----------------------------
# ggsave("seatac_tmax_1950_2021.png", p_max, width = 9, height = 5, dpi = 200)
# ggsave("seatac_tmin_1950_2021.png", p_min, width = 9, height = 5, dpi = 200)
# ggsave("seatac_tmean_1950_2021.png", p_mean, width = 9, height = 5, dpi = 200)
