#This code corresponds to the following publication
#Nitric and nitrous oxide fluxes from intensifying crop agriculture in the
#seasonally dry tropical Amazon-Cerrado border region
#Alexandra M. Huddell, Christopher Neill, Leonardo Maracahipes-Santos,
#Carlos Eduardo Pellegrino Cerri, and Duncan N. L. Menge

setwd("")

library(tidyverse)
library(ggplot2)
library(pracma)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(patchwork)
library(emmeans)


dat <- read.csv("NO_LU_flux_summary_11.19.19.csv")

# summarizing data by treatment, site, and date
dat <- as_tibble(dat)

#summarize data by date
treat_day_site_NO <- select(dat, day, month, year, site,
                            treatment, mean_flux) %>%
  filter(treatment != "harvest" & treatment != "planted") %>%
  group_by(treatment, site, day, month, year) %>%
  summarise(
    n = n(),
    NO_mg_N_m2_h = mean(mean_flux),
    sd = sd(mean_flux)
  )
treat_day_site_NO

treat_day_site_NO_all <- select(dat, day, month, year, site,
                                treatment, mean_flux) %>%
  filter(treatment != "harvest" & treatment != "planted")
treat_day_site_NO_all

#removing those levels from treatment and site
treat_day_site_NO$treatment <-
  as.character(treat_day_site_NO$treatment)
treat_day_site_NO$treatment <- as.factor(treat_day_site_NO$treatment)
levels(treat_day_site_NO$treatment)

treat_day_site_NO$site <- as.character(treat_day_site_NO$site)
treat_day_site_NO$site <- as.factor(treat_day_site_NO$site)
levels(treat_day_site_NO$site)

# convert date info in format 'mm/dd/yyyy'
treat_day_site_NO$date <-
  paste(
    as.character(treat_day_site_NO$month),
    as.character(treat_day_site_NO$day),
    as.character(treat_day_site_NO$year),
    sep = "/"
  )
treat_day_site_NO$date
treat_day_site_NO$date <-
  as.Date(treat_day_site_NO$date, "%m/%d/%Y")

#calculate standard error of daily fluxes
treat_day_site_NO$SE <- treat_day_site_NO$sd / sqrt(treat_day_site_NO$n)


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

#color palette
cbPalette <-
  c("#009E73", "#56B4E9", "#D55E00", "#F0E442", "#0072B2") #"#CC79A7","#E69F00"

#reorganize treatment levels
levels(treat_day_site_NO$treatment)
treat_day_site_NO$treatment <- ordered(treat_day_site_NO$treatment,
                                       levels = c("forest", "soy", "maize"))

#  facet label names
labs <- c("forest", "soybean", "soybean-maize")
names(labs) <- c("forest",  "soy", "maize")


LU_NO_plot <-
  ggplot(data = treat_day_site_NO, aes(x = date, y = NO_mg_N_m2_h)) +
  geom_jitter(
    aes(colour = treatment),
    width = .2,
    size = 3,
    alpha = .6
  ) +
  geom_line(aes(colour = factor(treatment))) +
  geom_errorbar(
    aes(
      ymax = NO_mg_N_m2_h + treat_day_site_NO$SE,
      ymin = NO_mg_N_m2_h - treat_day_site_NO$SE,
      colour = factor(treatment)
    )
  ) +
  ylab(expression(paste("NO flux (mg N ", " ", m ^ -2, " ", h ^ -1, ')', sep =
                          ""))) +
  xlab("Date") +
  theme_default() +
  geom_vline(xintercept = min(treat_day_site_NO$date) - 10) +
  geom_hline(yintercept = 0) +
  ylim(c(-0.25, 1)) +
  scale_color_manual(
    name = "Treatment",
    values = c(cbPalette[1], "blue", cbPalette[3]),
    breaks = c("forest", "soy", "maize"),
    labels = c("forest",  "soybean", "soybean-maize")
  ) +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  facet_grid( ~ treatment, labeller = labeller(treatment = labs)) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(breaks = as.Date(c(
    "2017-02-01", "2017-07-01",
    "2017-11-01",
    "2018-04-01"
  )),
  labels = c("Feb 17", "Jul 17", "Nov 17",
             "Apr 18"))

LU_NO_plot

#linear mixed effects model on NO flux by treatment and date
treat_day_site_NO$treatment = factor(treat_day_site_NO$treatment)
treat_day_site_NO$date = factor(treat_day_site_NO$date)
treat_day_site_NO$site = factor(treat_day_site_NO$site)

model = lmer(NO_mg_N_m2_h ~ treatment + date + (1 |
                                                  site), data = treat_day_site_NO)
anova(model)  #date is significant, not treatment
rand(model)   #random effects not significant
emmeans(model, specs = pairwise ~  "treatment")#post hoc test

# Calculating and plotting monthly fluxes -------------------------------------------------------------

#excluding harvest and planting treatments
dat <- filter(dat, treatment != "planted" & treatment != "harvest")
dat <- filter(dat, site != "planted" & treatment != "harvest")

#removing those levels from treatment and site
dat$treatment <- as.character(dat$treatment)
dat$treatment <- as.factor(dat$treatment)
levels(dat$treatment)

dat$site <- as.character(dat$site)
dat$site <- as.factor(dat$site)
levels(dat$site)

#converting dates to julian days then adds 365 to days in 2018 so they're chronological
# convert date info in format 'mm/dd/yyyy'
dat$date <- paste(as.character(dat$month),
                  as.character(dat$day),
                  as.character(dat$year),
                  sep = "/")
dat$date <- as.Date(dat$date, "%m/%d/%Y")

dat$j_day_raw <- as.numeric(format(as.Date(dat$date), "%j"))

as_tibble(dat)

#group by treatment, site, and day and then sample the fluxes randomly
chamb_lengths <- dat %>%
  group_by(site, treatment, j_day_raw) %>%
  mutate(unique_ID = paste0(site, treatment, j_day_raw)) %>%
  mutate(n_chamb = n()) #recording # of unique chambers per treatment/site/day combination

for (i in levels(as.factor(chamb_lengths$unique_ID))) {
  print (i)
}


sample_fluxes = NULL
# data.frame(site=factor(levels = c("T15", "T18", "T27")),
#                           treatment=factor(levels=c("80", "120", "160")),
#                           julian_day=integer(),
#                           replicate=factor(levels=c("a","b","c","d","e")),
#                           flux=numeric()
# )

#bootstrapping distributions of fluxes from each unique treatment, site, and sampling date to compile 5 flux observations for each sample event
set.seed(100)
flux_row_list <- list()
chamber_number <- 0
for (i in unique(chamb_lengths$unique_ID)) {
  chamber_number <- chamber_number + 1
  message(paste0("Chamber number:", chamber_number))
  tempdf = chamb_lengths[which(chamb_lengths$unique_ID == i), ]
  site = rep(tempdf$site[1], 5)
  treatment = rep(tempdf$treatment[1], 5)
  julian_day = rep(tempdf$j_day_raw[1], 5)
  replicate = c("a", "b", "c", "d", "e")
  
  
  if (tempdf$n_chamb[1] == 5) {
    print(paste("At", i, ": 5 chambers"))
    flux_pool = sample_n(as.data.frame(tempdf$mean_flux), 5, replace = F)
  } else if (tempdf$n_chamb[1] == 4) {
    print(paste("At", i, ": 4 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$mean_flux), 4, replace =
                             F),
                  sample_n(as.data.frame(tempdf$mean_flux), 1, replace = F))
  } else if (tempdf$n_chamb[1] == 3) {
    print(paste("At", i, ": 3 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$mean_flux), 3, replace =
                             F),
                  sample_n(as.data.frame(tempdf$mean_flux), 2, replace = F))
  } else if (tempdf$n_chamb[1] == 2) {
    print(paste("At", i, ": 2 chambers"))
    flux_pool = c(
      sample_n(as.data.frame(tempdf$mean_flux), 2, replace = F),
      sample_n(as.data.frame(tempdf$mean_flux), 2, replace = F),
      sample_n(as.data.frame(tempdf$mean_flux), 1, replace = F)
    )
  }
  
  #creating randomly sampled pools of 5 flux measurments with some repetition for those with n<5
  flux_pool = unname(unlist(flux_pool))
  
  #saving outputs
  flux_row <-
    data.frame(
      site = factor(
        x = site,
        levels = c(
          "app2",
          "area1",
          "area2",
          "mutum",
          "T15",
          "T18",
          "T27",
          "torre de soja"
        )
      ),
      treatment = factor(x = treatment, levels = c("forest", "maize", "soy")),
      julian_day = julian_day,
      replicate = replicate,
      flux = flux_pool
    )
  
  flux_row_list[[i]] <- flux_row
}

sample_fluxes <- as_tibble(bind_rows(flux_row_list))
sample_fluxes

#Now sampling each replicate/treatment/site combination and summing up cumulative flux
sample_fluxes$unique_ID <-
  paste0(sample_fluxes$site,
         sample_fluxes$treatment,
         sample_fluxes$replicate)
cumflux_row_list <- list()
for (i in unique(sample_fluxes$unique_ID)) {
  tempdf = filter(sample_fluxes, unique_ID == i) #subset each replicate by site and treatment
  tempdf <-
    tempdf[order(tempdf$julian_day), ] #order the julian days (x-axis for the trapezoidal integration in next step)
  cumflux = trapz(tempdf$julian_day, (tempdf$flux * 24))#multiply 24 h/d
  
  first_day = min(tempdf$julian_day)
  last_day = max(tempdf$julian_day)
  total_days = (last_day - first_day)
  
  #saving outputs
  cum_flux_row <-
    data.frame(
      site = factor(
        x = tempdf$site[1],
        levels = c(
          "app2",
          "area1",
          "area2",
          "mutum",
          "T15",
          "T18",
          "T27",
          "torre de soja"
        )
      ),
      treatment = factor(
        x = tempdf$treatment[1],
        levels = c("forest", "maize", "soy")
      ),
      replicate = tempdf$replicate[1],
      cumflux = cumflux,
      total_days = total_days[1]
    )
  
  cumflux_row_list[[i]] <- cum_flux_row
}


cumulative_fluxes <- as_tibble(bind_rows(cumflux_row_list))
cumulative_fluxes


cumulative_fluxes$kg_NO_N_ha <- cumulative_fluxes$cumflux / 10 ^ 2 *
  365 / cumulative_fluxes$total_days #converts from mg N/m2/day to kg N/ha/year

# Cumulative fluxes plots and stats ---------------------------------------

#values for table s2
cumulative_fluxes %>%
  group_by(treatment) %>%
  summarise(mean = mean(kg_NO_N_ha),
            SE = sd(kg_NO_N_ha) / sqrt(15))


cumulative_NO_fluxes <- cumulative_fluxes

#levels order
cumulative_NO_fluxes$site <- factor(
  cumulative_NO_fluxes$site,
  levels = c(
    "app2",
    "area1",
    "area2",
    "T15",
    "T18",
    "T27",
    "mutum",
    "torre de soja"
  )
)
cumulative_NO_fluxes$treatment <-
  ordered(cumulative_NO_fluxes$treatment,
          levels = c("forest",  "soy", "maize"))


annual_NO_plot <-
  ggplot(data = cumulative_NO_fluxes,
         aes(
           x = treatment,
           y = kg_NO_N_ha,
           color = site,
           shape = site,
           group = treatment
         ),) +
  geom_violin() +
  geom_jitter(size = 4, width = .05) +
  scale_shape_manual(name = "Site", values = (15:27)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.4) +
  #geom_segment(aes(x=0.1,xend=0.1,y=-5,yend=5), colour="black") +
  #geom_segment(aes(x=0,xend=3.5,y=0,yend=0),colour="black") +
  ylab(expression(paste("NO flux (kg NO-N  ",  ~ ha ^ -1, yr ^ -1, ")", sep = ""))) +
  xlab(" ") +
  scale_color_manual(
    name = "Site",
    values = RColorBrewer::brewer.pal(8, 'Dark2'),
    guide = guide_legend(ncol = 3)
  ) +
  theme_default() +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  annotate(
    "text",
    x = c(1, 2, 3),
    y = c(7, 7, 7),
    label = c("a", "b", "ab")
  ) +
  scale_x_discrete(labels = c("forest",  "soybean", "soybean-maize")) +
  theme(legend.background = element_rect(color = NA))


annual_NO_plot


#Combining NO and N2O cumulative fluxes
N2O_LU_cumulative_fluxes <- read.csv("N2O LU cumulative fluxes.csv")
N2O_LU_cumulative_fluxes$tracegas <-
  rep("N2O", nrow(N2O_LU_cumulative_fluxes))
cumulative_NO_fluxes$tracegas <- rep("NO", nrow(cumulative_NO_fluxes))
N2O_LU_cumulative_fluxes <-
  rename(N2O_LU_cumulative_fluxes, kg_N_ha = kg_N2O_N_ha) #changing column name
cumulative_NO_fluxes <-
  rename(cumulative_NO_fluxes, kg_N_ha = kg_NO_N_ha) #changing column name
LU_cumulative_fluxes <-
  rbind(N2O_LU_cumulative_fluxes, cumulative_NO_fluxes)
LU_cumulative_fluxes


cumulative_NO_fluxes %>% group_by(treatment) %>% summarise(mean = mean(kg_N_ha), SE =
                                                             sd(kg_N_ha) / sqrt(15))


#linear mixed effects model on cumulative fluxes ~ treatment
cumulative_NO_fluxes$treatment = factor(cumulative_NO_fluxes$treatment)

model2 = lmer(kg_N_ha ~ treatment + (1 |
                                       site), data = cumulative_NO_fluxes)
anova(model2)  #date is significant, not treatment
emmeans(model2, specs = pairwise ~  "treatment")#post hoc test

LU_cumulative_fluxes$treatment <-
  ordered(LU_cumulative_fluxes$treatment,
          levels = c("forest",  "soy", "maize"))
LU_cumulative_fluxes$tracegas <-
  ordered(LU_cumulative_fluxes$tracegas,
          levels = c("NO", "N2O"))

#plot with cumulative fluxes
annual_LU_plot <-
  ggplot(data = LU_cumulative_fluxes, aes(x = treatment, (y = kg_N_ha))) +
  geom_boxplot(aes(colour = factor(tracegas))) +
  ylab(expression(paste(
    N[2], "O or NO flux (kg N  ", ha ^ -1, yr ^ -1, ")", sep = ""
  ))) +
  xlab(" ") +
  theme_default() +
  geom_vline(xintercept = 0.4) +
  geom_hline(yintercept = 0) +
  scale_color_manual(
    name = expression("Trace gas"),
    labels = c(expression("NO", paste(N[2], "O"))),
    values = c("black", "red")
  ) +
  scale_x_discrete(labels = c("forest",  "soybean", "soybean-maize")) +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  annotate(
    "text",
    x = c(.8, 1.8, 2.8),
    y = c(6.2, 6.2, 6.2),
    label = c("a", "b", "ab")
  )
annual_LU_plot
