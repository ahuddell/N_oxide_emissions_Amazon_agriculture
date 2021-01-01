#This code corresponds to the following publication
#Nitric and nitrous oxide fluxes from intensifying crop agriculture in the
#seasonally dry tropical Amazon-Cerrado border region
#Alexandra M. Huddell, Christopher Neill, Leonardo Maracahipes-Santos,
#Carlos Eduardo Pellegrino Cerri, and Duncan N. L. Menge

setwd("")

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(lme4)
library(lmerTest)

#load data
fert_add_soil_N <-
  read.csv('fert_add_soil_N_data.csv',
           stringsAsFactors = F,
           header = T)

#transform to tidyverse object
soil_N <- as_tibble(fert_add_soil_N)


# summarizing by treatments and days post fertilization ----------------------------
#Sort by treatment and days post fertilization
fert_add_soil_N$treatment <- as.factor(fert_add_soil_N$treatment)
str(fert_add_soil_N$treatment)

#trimming dataset for required variables
NH4_mod <-
  select(fert_add_soil_N,
         days.post.fert,
         treatment,
         NH4_N_mg_kg_soil,
         site)
NO3_mod <-
  select(fert_add_soil_N,
         days.post.fert,
         treatment,
         NO3_N_mg_kg_soil,
         site)

#new code
NH4 <- NH4_mod %>% group_by(days.post.fert, treatment) %>%
  summarise(
    mean = mean(NH4_N_mg_kg_soil),
    n = n(),
    se = sd(NH4_N_mg_kg_soil) / sqrt(n())
  )

NO3 <- NO3_mod %>% group_by(days.post.fert, treatment) %>%
  summarise(
    mean = mean(NO3_N_mg_kg_soil),
    n = n(),
    se = sd(NO3_N_mg_kg_soil) / sqrt(n())
  )

# flux by day and treatment plot fert addition experiment------------------------------------------
theme_default <- function(axis_text_size = 13) {
  theme(
    text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = axis_text_size, colour = "black"),
    axis.text.y = element_text(size = axis_text_size, colour = "black"),
    axis.title.y = element_text(angle = 90, vjust = 0),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.key = element_blank(),
    legend.position = "bottom"
  )
}


#plot ammonium data
soil_fert_NH4 <-
  ggplot(NH4, aes(days.post.fert, mean, colour = treatment)) +
  geom_point(aes(color = treatment)) +
  geom_line(aes(colour = treatment)) +
  geom_errorbar(aes(
    ymin = (mean - se),
    ymax = (mean + se),
    width = .3
  )) +
  ylab(expression(paste(NH[4] ^ '+' ~ '(mg N kg' ~  ~ soil ^ {
    -1
  }, ')'))) +
  xlab('Time since fertilization (days)') +
  #ggtitle('Soil ammonium v. days since fertilization across N levels') +
  theme_default() +
  ylim(0, 200) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 200
  ), colour = "black") +
  geom_segment(aes(
    x = 0,
    xend = 30,
    y = 0,
    yend = 0
  ), colour = "black") +
  scale_color_manual(
    name = expression(paste("Treatment (kg N ", ha ^ -1, y ^ -1, ")   ", sep =
                              "")),
    values = c("gray55", "dodgerblue3", "midnight blue")
  ) +
  scale_shape_manual(name = expression(paste("Treatment (kg N ", ha ^ -1, y ^
                                               -1, ")   ", sep = "")),
                     values = c(16, 17, 18))
#plot nitrate data
soil_fert_NO3 <-
  ggplot(NO3, aes(days.post.fert, mean, colour = treatment)) +
  geom_point(aes(color = treatment)) +
  geom_line(aes(colour = treatment)) +
  geom_errorbar(aes(
    ymin = (mean - se),
    ymax = (mean + se),
    width = .3
  )) +
  ylab(expression(paste(NO[3] ^ '-' ~ '(mg N kg' ~  ~ soil ^ {
    -1
  }, ')'))) +
  xlab('Time since fertilization (days)') +
  #ggtitle('Soil nitrate v. days since fertilization across N levels') +
  theme_default() +
  ylim(0, 200) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 200
  ), colour = "black") +
  geom_segment(aes(
    x = 0,
    xend = 30,
    y = 0,
    yend = 0
  ), colour = "black") +
  scale_color_manual(
    name = expression(paste("Treatment (kg N ", ha ^ -1, y ^ -1, ")   ", sep =
                              "")),
    values = c("gray55", "dodgerblue3", "midnight blue")
  ) +
  scale_shape_manual(name = expression(paste("Treatment (kg N ", ha ^ -1, y ^
                                               -1, ")   ", sep = "")),
                     values = c(16, 17, 18))
soil_fert_NO3

#ANOVA on soil N with time

#reformatting predictors as factors
NH4_mod$treatment = factor(NH4_mod$treatment)
NH4_mod$days.post.fert = factor(NH4_mod$days.post.fert)
NH4_mod$site = factor(NH4_mod$site)
str(NH4_mod)

#linear mixed effect model on soil ammonium vs. days since fertilization and treatment
NH4_fert_mod <-
  (lmer(NH4_N_mg_kg_soil ~ days.post.fert + treatment + (1 | site),
        data = NH4_mod))
summary(NH4_fert_mod)
anova(NH4_fert_mod)

#forcing predictors to be factors
NO3_mod$treatment = factor(NO3_mod$treatment)
NO3_mod$days.post.fert = factor(NO3_mod$days.post.fert)
NO3_mod$site = factor(NO3_mod$site)
str(NO3_mod)

#linear mixed effect model on soil nitrate vs. days since fertilization and treatment
anova(lmer(NO3_N_mg_kg_soil ~ days.post.fert + treatment + (1 |
                                                              site),
           data = NO3_mod))
emmeans(NO3_fert_mod, specs = pairwise ~ "treatment")#post hoc test

# plots for LU data -------------------------------------------------------

#load land use data
LU_soil_N <- read.csv('LU_soil_N_data.csv', header = T)
LU_soil_N$date <- as.Date(LU_soil_N$date, format = "%m/%d/%Y")

#trimming dataset for required variables
NH4_LU_mod <- select(LU_soil_N, date, treatment, NH4_N_mg_kg_soil, site)
NO3_LU_mod <- select(LU_soil_N, date, treatment, NO3_N_mg_kg_soil, site)

#excluding harvest and planting treatments
NH4_LU_mod <-
  filter(NH4_LU_mod, treatment != "planted" & treatment != "harvest")
NO3_LU_mod <-
  filter(NO3_LU_mod, treatment != "planted" & treatment != "harvest")

#removing those levels from treatment and site
NH4_LU_mod$treatment <- as.character(NH4_LU_mod$treatment)
NH4_LU_mod$treatment <- as.factor(NH4_LU_mod$treatment)
levels(NH4_LU_mod$treatment)


NO3_LU_mod$treatment <- as.character(NO3_LU_mod$treatment)
NO3_LU_mod$treatment <- as.factor(NO3_LU_mod$treatment)
levels(NO3_LU_mod$treatment)


#summarize NH4 and NO3
NH4_LU <- NH4_LU_mod %>% group_by(date, treatment) %>%
  summarise(
    mean = mean(NH4_N_mg_kg_soil),
    n = n(),
    se = sd(NH4_N_mg_kg_soil) / sqrt(n())
  )

NO3_LU <- NO3_LU_mod %>% group_by(date, treatment) %>%
  summarise(
    mean = mean(NO3_N_mg_kg_soil),
    n = n(),
    se = sd(NO3_N_mg_kg_soil) / sqrt(n())
  )

#color palette
cbPalette <-
  c("#009E73", "#56B4E9", "#D55E00", "#F0E442", "#0072B2")

#reordering levels
levels(NH4_LU$treatment)
NH4_LU$treatment <-
  ordered(NH4_LU$treatment, levels = c("forest", "maize", "soy"))
levels(NO3_LU$treatment)
NO3_LU$treatment <-
  ordered(NO3_LU$treatment, levels = c("forest", "maize", "soy"))

#plot ammonium data
soil_LU_NH4 <- ggplot(NH4_LU, aes(date, mean, colour = treatment)) +
  geom_point(aes(color = treatment)) +
  geom_line(aes(colour = treatment)) +
  geom_errorbar(aes(
    ymin = (mean - se),
    ymax = (mean + se),
    width = .3
  )) +
  ylab(expression(paste(NH[4] ^ '+' ~ '(mg N kg' ~  ~ soil ^ {
    -1
  }, ')'))) +
  xlab('Date') +
  ylim(0, 60) +
  #ggtitle('Soil ammonium v. time across land uses') +
  theme_default() +
  geom_vline(xintercept = min(NH4_LU$date)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(
    name = "Treatment",
    values = c(cbPalette[1], cbPalette[3], "blue"),
    breaks = c("forest", "soy", "maize"),
    labels = c("forest",  "soybean", "soybean-maize")
  )


#plot nitrate data
soil_LU_NO3 <- ggplot(NO3_LU, aes(date, mean, colour = treatment)) +
  geom_point(aes(color = treatment)) +
  geom_line(aes(colour = treatment)) +
  geom_errorbar(aes(
    ymin = (mean - se),
    ymax = (mean + se),
    width = .3
  )) +
  ylab(expression(paste(NO[3] ^ '-' ~ '(mg N kg' ~  ~ soil ^ {
    -1
  }, ')'))) +
  xlab('Date') +
  #ggtitle('Soil nitrate v. time across land uses') +
  theme_default() +
  ylim(0, 60) +
  geom_vline(xintercept = min(NH4_LU$date)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(
    name = "Treatment",
    values = c(cbPalette[1], cbPalette[3], "blue"),
    breaks = c("forest", "soy", "maize"),
    labels = c("forest",  "soybean", "soybean-maize")
  )

#changing predictors to factors
NO3_LU_mod$treatment = factor(NO3_LU_mod$treatment)
NO3_LU_mod$date = factor(NO3_LU_mod$date)
NO3_LU_mod$site = factor(NO3_LU_mod$site)
str(NO3_LU_mod)

#linear mixed effect model on soil nitrate vs. date and treatment
NO3_LU_mod_fit <- (lmer(NO3_N_mg_kg_soil ~ date + treatment + (1 |
                                                                 site),
                        data = NO3_LU_mod))
summary(NO3_LU_mod_fit)
anova(NO3_LU_mod_fit)
emmeans(NO3_LU_mod_fit, specs = pairwise ~ "treatment")#post hoc test

#changing predictors to factors
NH4_LU_mod$treatment = factor(NH4_LU_mod$treatment)
NH4_LU_mod$days.post.fert = factor(NH4_LU_mod$date)
NH4_LU_mod$site = factor(NH4_LU_mod$site)
str(NH4_LU_mod)

#linear mixed effect model on soil ammonium vs. date and treatment
NH4_LU_mod_fit <- (lmer(NH4_N_mg_kg_soil ~ date + treatment + (1 |
                                                                 site),
                        data = NH4_LU_mod))
summary(NH4_LU_mod_fit)
anova(NH4_LU_mod_fit)
emmeans(NH4_LU_mod_fit, specs = pairwise ~ "treatment")#post hoc test
