# ===================================================================
#     MULTIPLE IMPUTATION (MICE) • COX GRAFT FAILURE • COX DEATH
# ===================================================================

# Packages --------------------------------------------------------------------
library(mice)
library(survival)
library(dplyr)
library(tidyr)
library(mitools)

# ===================================================================
# 1. LOAD DATA + PREPARE FOR IMPUTATION
# ===================================================================

df_original <- read.csv(
  "C:/Users/tyenc/OneDrive - UvA/Documenten/SCHOOL/UVA/Courses/3. Advanced Data Analysis in Medicine/Week 4/d.csv"
)

# Select relevant variables for imputation
df_impute <- df_original %>%
  select(
    ID, age_at_tx, sex_pat, bmi, sexdon, agedon, type_dia, duur_dia, retrans,
    gfr1:gfr7, screat1:screat7, map1:map7,
    stat_gra, stat_pat, time_to_graft_failure, time_to_death
  ) %>%
  mutate(
    sex_pat  = factor(sex_pat),
    sexdon   = factor(sexdon),
    type_dia = factor(type_dia),
    stat_gra = factor(stat_gra),
    stat_pat = factor(stat_pat)
  )

# ===================================================================
# 2. RESTRICTED MICE IMPUTATION (ONLY BIOMARKERS)
# ===================================================================

ini  <- mice(df_impute, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# Identify biomarkers
lab_cols <- names(df_impute)[grepl("gfr|screat|map", names(df_impute))]

# Variables NOT imputed
no_impute <- setdiff(names(df_impute), lab_cols)
meth[no_impute] <- ""
pred[, no_impute] <- 0
pred[no_impute, ] <- 0

# Biomarker imputation method
meth[lab_cols] <- "pmm"
pred[lab_cols, lab_cols] <- 1

set.seed(42)
imp <- mice(df_impute, m = 10, maxit = 10, method = meth, predictorMatrix = pred)

# ===================================================================
# 3. FUNCTION: Wide → Long (Biomarkers)
# ===================================================================

make_long <- function(df) {
  df %>%
    pivot_longer(
      cols = c(gfr1:gfr7, screat1:screat7, map1:map7),
      names_to = c("marker", "tp"),
      names_pattern = "([a-z]+)([1-7])",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = marker,
      values_from = value
    ) %>%
    mutate(
      tp = as.integer(tp),
      years = case_when(
        tp == 1 ~ 0.25,
        tp == 2 ~ 0.50,
        tp == 3 ~ 1.00,
        tp == 4 ~ 2.00,
        tp == 5 ~ 5.00,
        tp == 6 ~ 10.00,
        tp == 7 ~ 15.00
      )
    )
}

# Convert all imputations to long format
all_long <- lapply(1:imp$m, function(i) make_long(complete(imp, i)))

# ===================================================================
# 4. PROCESS LONG DATA FOR TIME-DEPENDENT COX (GRAFT FAILURE)
# ===================================================================

process_td <- function(df) {
  df %>%
    arrange(ID, years) %>%
    group_by(ID) %>%
    mutate(
      censor_time = pmin(time_to_graft_failure, time_to_death),
      tstart_raw  = lag(years, default = 0),
      tstop_raw   = years,
      tstart      = pmin(tstart_raw, censor_time),
      tstop       = pmin(tstop_raw, censor_time)
    ) %>%
    filter(tstart < censor_time) %>%
    mutate(
      graft_event = ifelse(
        stat_gra == "Graft loss" &
          time_to_graft_failure > tstart &
          time_to_graft_failure <= tstop,
        1, 0
      )
    ) %>%
    ungroup()
}

dlong_td_list <- lapply(all_long, process_td)

# ===================================================================
# 5. FIT COX MODELS Baseline Cox model (no time-dependent covariates)
# ===================================================================
library(survival)
library(mitools)

# Event indicator
df_impute$event_gra <- as.numeric(df_impute$stat_gra == "Graft loss")

# Fit baseline Cox model on each imputed dataset
cox_baseline_fits <- lapply(1:imp$m, function(i) {
  df <- complete(imp, i)
  df$event_gra <- as.numeric(df$stat_gra == "Graft loss")
  
  coxph(
    Surv(time_to_graft_failure, event_gra) ~
      age_at_tx + sex_pat + bmi + type_dia + duur_dia +
      sexdon + agedon + retrans,
    data = df
  )
})

# ===================================================================
# RUBIN'S POOLING — GRAFT FAILURE
# ===================================================================


# Pool results
pooled_baseline <- MIcombine(cox_baseline_fits)
pooled_baseline

b <- pooled_baseline$coefficients
v <- diag(pooled_baseline$variance)

pooled_table_baseline <- data.frame(
  Variable = names(b),
  HR       = exp(b),
  Lower95  = exp(b - 1.96*sqrt(v)),
  Upper95  = exp(b + 1.96*sqrt(v)),
  pvalue   = 2*(1 - pnorm(abs(b/sqrt(v))))
)

pooled_table_baseline

# ===================================================================
# 6. Cox model with time-dependent covariates
# ===================================================================

cox_tdc_fits <- lapply(1:imp$m, function(i) {
  df_wide <- complete(imp, i)
  df_long <- make_long(df_wide)
  
  df_wide$event_gra <- as.numeric(df_wide$stat_gra == "Graft loss")
  df_wide$gfstop    <- df_wide$time_to_graft_failure
  
  df_tdc <- tmerge(
    data1 = df_wide,
    data2 = df_wide,
    id = ID,
    tstop = gfstop,
    graft_event = event(time_to_graft_failure, event_gra)
  )
  
  df_tdc <- tmerge(
    data1 = df_tdc,
    data2 = df_long,
    id = ID,
    gfr_tdc    = tdc(years, gfr),
    screat_tdc = tdc(years, screat),
    map_tdc    = tdc(years, map)
  )
  
  coxph(
    Surv(tstart, tstop, graft_event) ~
      age_at_tx + sex_pat + bmi + type_dia + duur_dia +
      sexdon + agedon + retrans +
      gfr_tdc + screat_tdc + map_tdc,
    data = df_tdc,
    id = ID
  )
})

# ===================================================================
# RUBIN'S POOLING — GRAFT FAILURE
# ===================================================================

pooled_tdc <- MIcombine(cox_tdc_fits)

b <- pooled_tdc$coefficients
v <- diag(pooled_tdc$variance)

pooled_table_tdc <- data.frame(
  Variable = names(b),
  HR       = exp(b),
  Lower95  = exp(b - 1.96*sqrt(v)),
  Upper95  = exp(b + 1.96*sqrt(v)),
  pvalue   = 2*(1 - pnorm(abs(b/sqrt(v))))
)

pooled_table_tdc

# ===================================================================
# plot
# ===================================================================
# Prepare data
forest_tdc <- pooled_table_tdc %>%
  mutate(
    Variable = factor(Variable, levels = rev(Variable))  # reverse order for plot
  )

# Forest plot
ggplot(forest_tdc, aes(x = HR, y = Variable)) +
  geom_point(size = 3, color = "#1F618D") +
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.2, color = "#1F618D") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
  scale_x_log10() +
  theme_minimal(base_size = 13) +
  xlab("Hazard Ratio (log scale)") +
  ylab("") +
  ggtitle("Forest Plot: Cox Model with Time-Dependent Covariates for Graft Failure")
# replace scale_x_log10() with:
  scale_x_continuous()
  