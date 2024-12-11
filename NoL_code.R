### TMS Aphasia PNT Analysis - 6.20.2024 ###

#----Set-up----
library(dplyr); library(ppcor); library(emmeans); library(sjPlot); library(ggplot2); library(lme4); library(lmerTest); library(sjPlot)
pnt_df <- read.csv("TMSA_PNT_n30_100124.csv") 
pnt_df <- pnt_df %>% mutate_at(c(5,11:16), as.numeric)
pnt_df <- pnt_df %>%mutate_at(c(1:4,6,8,10), factor)

#----Data cleaning----
pnt_df <- pnt_df[!(pnt_df$TrialDrop=="1"),]
pnt_df <- pnt_df[!(pnt_df$SubjID=="TMSA-04"),] #unintelligible
pnt_df <- pnt_df[!(pnt_df$SubjID=="TMSA-24"),] #"has-to" perseveration

#----Contrast setting----
contrasts(pnt_df$Stim) <- c(1,0)
contrasts(pnt_df$TimePoint) <- c(1,0)
contrasts(pnt_df$TimePoint) <- factor(c("Baseline", "3mo", "6mo"))
contrasts(time) <- matrix(c(-1, 1, 0,  # baseline vs 3-months
                            -1, 0, 1   # baseline vs 6-months
                          ), ncol = 2)
colnames(attr(pnt_df$Stim, "contrasts")) <- c("TMS")

#----LM models-----
lmSwt <- lm(ProporImprv ~ 1 + Stim*TimePoint + Stim*BaselineSwt, data=pnt_df)
lmPwt <- lm(ProporImprv ~ 1 + Stim*TimePoint + Stim*BaselinePwt, data=pnt_df)

# effect sizes
if (!require(effectsize)) {
  install.packages("effectsize")
}
library(effectsize)
eta_squared(lmSwt, partial = TRUE)
eta_squared(lmPwt, partial = TRUE)

# VIF
vif(lmSwt)
vif(lmPwt)
lmBothWt <- lm(ProporImprv ~ 1 + Stim*TimePoint + Stim*BaselineSwt + Stim*BaselinePwt, data=pnt_df)
vif(lmBothWt)

#----Severity correlations-----
cor.test(pnt_df$BaselineWABAQ, pnt_df$BaselineSwt, method="pearson") # r=0.66, p<0.001
cor.test(pnt_df$BaselineWABAQ, pnt_df$BaselinePwt, method="pearson") # r=0.34, p<0.001
cor.test(pnt_df$BaselineSwt, pnt_df$BaselinePwt, method="pearson") # r=0.11, p<0.001

#----LOOCV for over-fitting-----
train_control_loocv <- trainControl(method = "LOOCV")
Scv_model_loocv <- train(
  Smodel_formula,
  data = pnt_df_clean,
  method = "lm",
  trControl = train_control_loocv)
print(Scv_model_loocv)
Pcv_model_loocv <- train(
  Pmodel_formula,
  data = pnt_df_clean,
  method = "lm",
  trControl = train_control_loocv)
print(Pcv_model_loocv)

#----absolute outcomes (appendix)-----
pnt_df$TimePoint <- factor(pnt_df$TimePoint, levels = c("Baseline", "3mo", "6mo"))
contrasts(pnt_df$TimePoint) <- matrix(c(-1, 1, 0,  # Baseline vs 3mo
                                        -1, 0, 1), # Baseline vs 6mo
                                        ncol = 2)
colnames(contrasts(pnt_df$TimePoint)) <- c("Baseline_vs_3mo", "Baseline_vs_6mo")

lmSwt_abs <- lm(Acc_Bin ~ 1 + Stim*TimePoint + Stim*BaselineSwt, data=pnt_df)
lmPwt_abs <- lm(Acc_Bin ~ 1 + Stim*TimePoint + Stim*BaselinePwt, data=pnt_df)

eta_squared(lmSwt_abs, partial = TRUE)
eta_squared(lmPwt_abs, partial = TRUE)

emmSwt <- emmeans(lmSwt_abs, ~ Stim*TimePoint + Stim*BaselineSwt)
emmPwt <- emmeans(lmPwt_abs, ~ Stim*TimePoint + Stim*BaselinePwt)
