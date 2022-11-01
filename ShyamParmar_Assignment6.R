# Clear Console
rm(list=ls())

# install.packages("pacman")
pacman::p_load(dplyr, tidyr, caret, ggplot2, caTools, corrplot, PerformanceAnalytics, AER, MASS, stargazer, pscl, jtools, Hmisc, ggcorrplot, rpart, rpart.plot, survival, survminer)

# Read file
df <- read.table('C:/Users/Scott/Downloads/survival.txt', skip = 15)

# Rename the headers
names(df) <- c('treatment', 'cell_type', 'survival', 'status', 'ksrnofsky_score', 'months_from_diag', 'age', 'prior_chemo')

# Convert prior_chemo to binary
df$prior_chemo <- ifelse(df$prior_chemo == 0, 0, 1)

# Descriptive analysis
unique(df$cell_type)
table(df$status)
table(df$prior_chemo)
hist(df$survival)
hist(df$ksrnofsky_score)
hist(df$months_from_diag)
hist(df$age)
table(df$treament)
chart.Correlation(df[as.integer(which(sapply(df,class)=="integer"))]) # Plot for integer variables

# Create age group bins
df <- df %>% mutate(agegroup = case_when(
  age >= 80  & age <= 89 ~ '6',
  age >= 70  & age <= 79 ~ '5',
  age >= 60  & age <= 69 ~ '4',
  age >= 50  & age <= 59 ~ '3',
  age >= 40  & age <= 49 ~ '2',
  age >= 30  & age <= 39 ~ '1'))

df$treatment <- as.factor(df$treatment)

# K 
y <- Surv(df$months_from_diag, df$status)              # Y is a combination of time and event

km1 <- survfit(y ~ 1)
summary(km1)
plot(km1, xlab="Time", ylab="Survival Probability")

# Kaplan-Meier non-parametric analysis by group
km2 <- survfit(y ~ df$treatment)
summary(km2)
plot(km2, xlab="Time", ylab="Survival Probability")

fit <- survfit(Surv(df$months_from_diag, df$status) ~ treatment, data = df)
ggsurvplot(fit, xlab="Time (Month)", break.x.by = 5, conf.int = TRUE, surv.median.line = "hv", xlim = c(0, 45), data = df)

#############################################################################

# Cox proportional hazard model - coefficients and hazard rates
cox <- coxph(y ~ df$treatment + df$agegroup + df$prior_chemo)
summary(cox)

# Exponential, Weibull, and log-logistic parametric model coefficients
exp <- survreg(y ~ df$treatment + df$agegroup + df$prior_chemo, dist="exponential")
summary(exp)

weibull <- survreg(y ~ df$treatment + df$agegroup + df$prior_chemo, dist="weibull")
summary(weibull)

loglogistic <- survreg(y ~ df$treatment + df$age + df$prior_chemo, dist="loglogistic")
summary(loglogistic)

library(stargazer)
stargazer(cox, exp, weibull, type="text", single.row=TRUE)


