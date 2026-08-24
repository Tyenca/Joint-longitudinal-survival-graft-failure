# ===================================================================
#   LME MODEL COMPARISON + JOINT MODELS (GFR / Creat / MAP)
# ===================================================================

library(nlme)
library(splines)
library(survival)
library(JM)

d_long <- dlong
d_base <- d
# ===================================================================
# 1. FUNCTION: Fit 3 LME Models (RI / RS / Spline Random)
# ===================================================================

fit_lme_models <- function(df, biomarker) {
  
  # Fixed effects: spline in years
  formula_fixed <- as.formula(
    paste0(biomarker, " ~ ns(years, df = 3)")
  )
  
  # 1. Random intercept
  m1 <- lme(
    formula_fixed,
    random = ~ 1 | ID,
    data = df,
    na.action = na.exclude,
    control = lmeControl(msMaxIter = 200)
  )
  
  # 2. Random intercept + random slope in time
  m2 <- lme(
    formula_fixed,
    random = ~ years | ID,
    data = df,
    na.action = na.exclude,
    control = lmeControl(msMaxIter = 200)
  )
  
  # 3. Random spline (df = 2)
  m3 <- lme(
    formula_fixed,
    random = ~ ns(years, df = 2) | ID,
    data = df,
    na.action = na.exclude,
    control = lmeControl(msMaxIter = 200)
  )
  
  # Compare AICs
  aic_table <- AIC(m1, m2, m3)
  
  list(
    models = list(m1 = m1, m2 = m2, m3 = m3),
    AIC = aic_table
  )
}

# ===================================================================
# 2. LME MODEL SELECTION FOR EACH BIOMARKER
# ===================================================================

lme_gfr_res   <- fit_lme_models(dlong, "gfr")
lme_creat_res <- fit_lme_models(dlong, "creat")
lme_map_res   <- fit_lme_models(dlong, "map")

lme_gfr_res$AIC
lme_creat_res$AIC
lme_map_res$AIC


# ===================================================================
# 3. PREPARE LONGITUDINAL DATASETS
# ===================================================================

d_long_gfr   <- d_long[, c("ID", "years", "gfr")]
d_long_creat <- d_long[, c("ID", "years", "creat")]
d_long_map   <- d_long[, c("ID", "years", "map")]

# ===================================================================
# 4. FIT FINAL LME MODELS (Random Intercept + Spline Fixed)
# ===================================================================

lme_gfr <- lme(
  gfr ~ ns(years, df = 3),
  random = ~ 1 | ID,
  data = d_long_gfr,
  na.action = na.exclude,
  control = lmeControl(msMaxIter = 200)
)


lme_creat <- lme(
  creat ~ ns(years, df = 3),
  random = ~ ns(years, df = 2) | ID,   # ⬅️ randomized spline terms
  data = d_long_creat,
  na.action = na.exclude,
  control = lmeControl(msMaxIter = 200)
)


lme_map <- lme(
  map ~ ns(years, df = 3),
  random = ~ ns(years, df = 2) | ID,
  data = d_long_map,
  na.action = na.exclude,
  control = lmeControl(msMaxIter = 200)
)


# ===================================================================
# 5. COX MODEL FOR SURVIVAL PROCESS
# ===================================================================

cox_graft <- coxph(
  Surv(time_to_graft_failure, stat_gra == "Graft loss") ~ 
    age_at_tx + sex_pat + bmi + type_dia + duur_dia +
    sexdon + agedon + retrans,
  data  = d_base,     # baseline dataset
  x     = TRUE,
  model = TRUE
)

# ===================================================================
# 6. FIT JOINT MODELS (Weibull PH)
# ===================================================================

joint_gfr <- jointModel(
  lmeObject  = lme_gfr,
  survObject = cox_graft,
  timeVar    = "years",
  method     = "weibull-PH-aGH"
)

joint_creat <- jointModel(
  lmeObject  = lme_creat,
  survObject = cox_graft,
  timeVar    = "years",
  method     = "weibull-PH-aGH"
)

joint_map <- jointModel(
  lmeObject  = lme_map,
  survObject = cox_graft,
  timeVar    = "years",
  method     = "weibull-PH-aGH"
)

summary(joint_gfr)
summary(joint_creat)
summary(joint_map)

# ===================================================================
# 7. Longitudinal Trajectory Plot (Predicted GFR / creat / MAP over time)
# ===================================================================

plot(joint_gfr, param = "Longitudinal")
plot(joint_creat, param = "Longitudinal")
plot(joint_map, param = "Longitudinal")

# ===================================================================
# 8. Predicted Hazard Over Time
# ===================================================================

plot(joint_gfr, param = "Event")
plot(joint_creat, param = "Event")
plot(joint_map, param = "Event")

