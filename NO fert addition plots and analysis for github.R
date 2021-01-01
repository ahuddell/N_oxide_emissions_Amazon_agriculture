#This code corresponds to the following publication
#Nitric and nitrous oxide fluxes from intensifying crop agriculture in the
#seasonally dry tropical Amazon-Cerrado border region
#Alexandra M. Huddell, Christopher Neill, Leonardo Maracahipes-Santos,
#Carlos Eduardo Pellegrino Cerri, and Duncan N. L. Menge

setwd("")

dat <- read.csv("fert_add_flux_summary_7.21.20.csv")
precip <- read.csv("Tanguro_precip.csv")

library(tidyverse)
library(ggplot2)
library(pracma)
library(lme4)
library(lmerTest)
library(emmeans)

# Plot data ---------------------------------------------------------------

# summarizing data by treatment, site, and day post fertilizer
dat <- as_tibble(dat)

treat_day_site_NO <- dat %>%
  select(days.post.fert, site, treatment, mean_flux, lowb, upb) %>%
  group_by(treatment, site, days.post.fert) %>%
  summarise(
    n = n(),
    NO_mg_N_m2_h = mean(mean_flux),
    sd = sd(mean_flux),
    lowb_indiv_flux = min(lowb)
  )


#calculating standard errors on daily NO data
treat_day_site_NO$SE <-
  treat_day_site_NO$sd / sqrt(treat_day_site_NO$n)

treat_day_site_NO_all <- select(dat, days.post.fert, site,
                                treatment, mean_flux) %>%
  filter(treatment != "harvest" & treatment != "planted")
treat_day_site_NO_all


# flux by day and treatment plot ------------------------------------------
#theme for plots
theme_default <- function(axis_text_size = 13) {
  theme(
    text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = axis_text_size, colour = "black"),
    axis.text.y = element_text(size = axis_text_size, colour = "black"),
    axis.title.y = element_text(angle = 90, vjust = 0),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.position = "none"
  )
}

#ordering sites
levels(treat_day_site_NO$site)
treat_day_site_NO$site <- ordered(treat_day_site_NO$site,
                                  levels = c("T15", "T18", "T27"))

#labels for facet
labs <- c("site 1", "site 2", "site 3")
names(labs) <- c("T15", "T18", "T27")


NO_fert <-
  ggplot() +
  geom_jitter(
    data = treat_day_site_NO_all,
    aes(x = days.post.fert, y = mean_flux),
    size = 3,
    width = 2,
    color = "darkgrey",
    alpha = .8
  ) +
  geom_point(
    data = treat_day_site_NO,
    aes(
      x = days.post.fert,
      y = NO_mg_N_m2_h,
      colour = factor(treatment)
    ),
    shape = 17,
    size = 4
  ) +
  geom_line(data = treat_day_site_NO, aes(
    x = days.post.fert,
    y = NO_mg_N_m2_h,
    colour = factor(treatment)
  )) +
  geom_errorbar(
    data = treat_day_site_NO,
    aes(
      x = days.post.fert,
      ymax = NO_mg_N_m2_h + treat_day_site_NO$SE,
      ymin = NO_mg_N_m2_h - treat_day_site_NO$SE,
      colour = factor(treatment)
    )
  ) +
  ylab(expression(paste("NO flux (mg N ", " ", m^-2, " ", h^-1, ")",
    sep =
      ""
  ))) +
  xlab("Time since fertilization (days)") +
  theme_default() +
  geom_segment(
    data = treat_day_site_NO, 
    aes(x = 0, xend = 0,
        y = 0, yend = 6), 
    colour = "black"
  ) +
  geom_segment(
    data = treat_day_site_NO, 
    aes(
      x = 0, xend = 32,
      y = 0, yend = 0
    ), 
    colour = "black"
  ) +
  scale_color_manual(
    name = expression(paste("Treatment (kg N ", ha^-1, yr^
      -1, ")   ", sep = "")),
    values = c("deepskyblue2", "dodgerblue3", "navy")
  ) +
  scale_shape_manual(
    name = expression(paste("Treatment (kg N ", ha^-1, yr^
      -1, ")   ", sep = "")),
    values = c(16, 17, 18)
  ) +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  facet_grid(~site, labeller = labeller(site = labs)) +
  theme(strip.background = element_rect(fill = "white"))
NO_fert

#mixed effects model on daily NO fluxes by treatment and days since fertilization
treat_day_site_NO$treatment = factor(treat_day_site_NO$treatment)
treat_day_site_NO$days.post.fert = factor(treat_day_site_NO$days.post.fert)
treat_day_site_NO$site = factor(treat_day_site_NO$site)

model = lmer(NO_mg_N_m2_h ~ treatment + days.post.fert + (1 | site), 
             data = treat_day_site_NO)
anova(model)  #date and treatment significant
ranef(model)
lmerTest::rand(model)   #random effects not significant
emmeans(model, specs = pairwise ~  "treatment")#post hoc test

#precipitation plot
precip %>%
  filter(days.post.fert < 8) %>%
  group_by(site) %>%
  summarize(sum = round(sum(precip_mm), 0))

precip_names <- list("T15" = "site 1\n\ndays 0-7 total 64mm",
                     "T18" = "site 2\n\ndays 0-7 total 53mm",
                     "T27" = "site 3\n\ndays 0-7 total 130mm")

precip_labeller <- function(variable, value) {
  return(precip_names[value])
}

precip_plot <- ggplot(data = precip,
                      aes(x = days.post.fert,
                          y = precip_mm)) +
  geom_col(alpha = .8) +
  ylab("Precipitation (mm)") +
  xlab("Time since fertilzation (days) ") +
  theme_default() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

precip_plot <-
  precip_plot + facet_grid(. ~ site, labeller = precip_labeller) +
  theme(strip.background = element_rect(fill = "white"))

precip_plot


# Calculating and plotting monthly fluxes -------------------------------------------------------------

#group by treatment, site, and day and then sample the fluxes randomly
chamb_lengths <- dat %>%
  group_by(site, treatment, days.post.fert) %>%
  mutate(unique_ID = paste0(site, treatment, days.post.fert)) %>%
  mutate(n_chamb = n()) #recording # of unique chambers per treatment/site/day combination


sample_fluxes = data.frame(
  site = factor(levels = c("T15", "T18", "T27")),
  treatment = factor(levels = c("80", "120", "160")),
  days.post.fert = integer(),
  replicate = factor(levels = c("a", "b", "c", "d", "e")),
  flux = numeric()
)
set.seed(100)
flux_row_list <- list()
chamber_number <- 0


#bootstrapping distributions of fluxes from each unique treatment, site, and sampling date to compile 5 flux observations for each sample event
for (i in unique(chamb_lengths$unique_ID)) {
  chamber_number <- chamber_number + 1
  message(paste0("Chamber number:", chamber_number))
  tempdf = chamb_lengths[which(chamb_lengths$unique_ID == i),]
  site = rep(tempdf$site[1], 5)
  treatment = rep(tempdf$treatment[1], 5)
  days.post.fert = rep(tempdf$days.post.fert[1], 5)
  replicate = c("a", "b", "c", "d", "e")
  
  if (tempdf$n_chamb[1] == 5) {
    print(paste("At", i, ": 5 chambers"))
    flux_pool = sample_n(as.data.frame(tempdf$mean_flux), 5, replace =
                           F)
  } else if (tempdf$n_chamb[1] == 4) {
    print(paste("At", i, ": 4 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$mean_flux), 4, replace =
                             F),
                  sample_n(as.data.frame(tempdf$mean_flux), 1, replace =
                             F))
  } else if (tempdf$n_chamb[1] == 3) {
    print(paste("At", i, ": 3 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$mean_flux), 3, replace =
                             F),
                  sample_n(as.data.frame(tempdf$mean_flux), 2, replace =
                             F))
  }
  
  #creating randomly sampled pools of 5 flux measurments with some repetition for those with n<5
  flux_pool = unname(unlist(flux_pool))
  
  #saving outputs
  flux_row <-
    data.frame(
      site = factor(x = site, levels = c("T15", "T18", "T27")),
      treatment = factor(x = treatment, levels = c("80", "120", "160")),
      days.post.fert = days.post.fert,
      replicate = replicate,
      flux = flux_pool
    )
  
  flux_row_list[[i]] <- flux_row
}

sample_fluxes <- as_tibble(bind_rows(flux_row_list))

#Now sampling each replicate/treatment/site combination and summing up cumulative flux
sample_fluxes$unique_ID <-
  paste0(sample_fluxes$site,
         sample_fluxes$treatment,
         sample_fluxes$replicate)
cumflux_row_list <- list()

for (i in unique(sample_fluxes$unique_ID)) {
  tempdf = filter(sample_fluxes, unique_ID == i) #subset each replicate by site and treatment
  tempdf <- tempdf[order(tempdf$days.post.fert),]
  cumflux = trapz(tempdf$days.post.fert, (tempdf$flux * 24))#multiply 24 h/d
  max_day = max(tempdf$days.post.fert)
  
  #saving outputs
  cum_flux_row <-
    data.frame(
      site = factor(x = tempdf$site[1], levels = c("T15", "T18", "T27")),
      treatment = factor(
        x = tempdf$treatment[1],
        levels = c("80", "120", "160")
      ),
      days.post.fert = tempdf$days.post.fert[1],
      replicate = tempdf$replicate[1],
      cumflux = cumflux,
      max_day = max_day
    )
  
  cumflux_row_list[[i]] <- cum_flux_row
}

cumulative_fluxes <- as_tibble(bind_rows(cumflux_row_list))
cumulative_fluxes

#calculating cumulative monthly fluxes 
cumulative_fluxes$kg_NO_N_ha <- cumulative_fluxes$cumflux / 10 ^ 2 *
  30 / cumulative_fluxes$max_day #converts from mg N/m2/day to kg N/ha/month

write.csv(cumulative_fluxes,"NO_fertadd_cumulative_flux_month.csv",row.names = F)

#mixed effects model on fluxes ~ treatment
cumulative_fluxes$treatment = factor(cumulative_fluxes$treatment)

model2 = lmer(kg_NO_N_ha ~ treatment + (1 |
                                          site), data = cumulative_fluxes)
anova(model2)  #date is significant, not treatment
emmeans(model2, specs = pairwise ~  "treatment")#post hoc test



#values for table s1
cumulative_fluxes %>%
  group_by(treatment) %>%
  summarize(mean = mean(kg_NO_N_ha),
            SE = sd(kg_NO_N_ha) / sqrt(5))


monthly_NO_plot <- ggplot(data = cumulative_fluxes,
                          aes(x = treatment,
                              y = kg_NO_N_ha,
                              color = treatment), ) +
  geom_jitter(size = 4,
              width = .05,
              alpha = .8) +
  scale_shape_manual(name = "Site", values = (15:27)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.4) +
  ylab(expression(paste("NO flux (kg N ", ha ^ -1, month ^ -1, ")", sep = ""))) +
  xlab(expression(paste("Treatment (kg N ", ha ^ -1, yr ^ -1, ")", sep =
                          ""))) +
  scale_color_manual(
    name = expression(paste("Treatment (kg N ", ha ^ -1, yr ^ -1, ")   ", sep =
                              "")),
    values = c("gray55", "dodgerblue3", "midnight blue")
  ) +
  theme_default() +
  theme(legend.position = "bottom", legend.key = element_blank())
NO_fert_cumulative <-
  monthly_NO_plot + facet_grid(~ site) + theme(strip.background = element_rect(fill =
                                                                                 "white"))
NO_fert_cumulative



# Cumulative NO and N2O plot--------------------------------------------

#reading in N2O data and combining with cumulative NO fluxes 
n2o <- read.csv("N2O_fertadd_cumulative_flux_month.csv")

cumul_NO <-
  select(cumulative_fluxes, site, treatment, kg_N_ha = kg_NO_N_ha)
cumul_NO$trace_gas <- rep("NO", 45)

cumul_N2O <- select(n2o, site, treatment, kg_N_ha = kg_N2O_N_ha)
cumul_N2O$trace_gas <- rep("N2O",45)

N2O_NO <- as_tibble(rbind(cumul_N2O, cumul_NO))
N2O_NO$treatment <- as.numeric(N2O_NO$treatment)

#summarize cumulative 
N2O_NO_summary <- N2O_NO %>%
  group_by(trace_gas, treatment) %>%
  summarize(mean = mean(kg_N_ha), SE = (sd(kg_N_ha) / sqrt(15)))

#reordering trace gases
N2O_NO_summary$trace_gas <-
  as.factor(ordered(N2O_NO_summary$trace_gas, c("NO", "N2O")))

#plot
NO_N2O_cumulative <-
  ggplot(data = N2O_NO_summary,
         aes(
           x = as.numeric(treatment),
           y = mean,
           ymax = mean + SE,
           ymin = mean - SE,
           colour = factor(trace_gas)
         )) +
  geom_pointrange(position = position_jitter(width = .1), alpha = .7) +
  geom_line(aes(colour = factor(trace_gas)), alpha = .6) +
  ylab(expression(paste(
    N[2], "O or NO flux (kg N  ", ha ^ -1, mo ^ -1, ")", sep = ""
  ))) +
  xlab(expression(paste("Treatment (kg N ", ha ^ -1, yr ^ -1, ")", sep =
                          ""))) +
  theme_default() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 75) +
  scale_color_manual(
    name = expression("Trace gas"),
    labels = c(expression("NO", paste(N[2], "O"))),
    values = c("black", "red")
  ) +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    strip.background = element_rect(fill = "white")
  ) +
  annotate(
    "text",
    x = c(80, 120, 160),
    y = c(8.2, 8.2, 8.2),
    label = c("a", "ab", "c")
  ) +
  annotate(
    "text",
    x = c(80, 120, 160),
    y = c(4, 4, 4),
    label = c("a", "b", "c"),
    col = "red"
  )
NO_N2O_cumulative

