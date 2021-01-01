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
library(emmeans)

dat <- read.csv("N2O_flux_summary_11.12.19.csv")
precip <- read.csv("Rainfall_Estacao_Darro_2013_2019.csv")


## summarizing data by treatment, site, and day post fertilizer
dat <- tbl_df(dat)
treat_day_site_N2O <-
  dat %>%
  filter(LU == 1) %>%
  filter(treatment != "harvest" & treatment != "planted") %>%
  select(day, month, year, site, treatment, chamber, N2O_mg_N_m2_h) %>%
  group_by(treatment, site, day, month, year) %>%
  summarise(
    n = n(),
    mean_N2O_mg_N_m2_h = mean(N2O_mg_N_m2_h, na.rm = T),
    sd = sd(N2O_mg_N_m2_h, na.rm = T)
  )

treat_day_site_N2O_all <-
  dat %>%
  filter(LU == 1) %>%
  filter(treatment != "harvest" & treatment != "planted")

#organizing treatment levels
treat_day_site_N2O$treatment <-
  as.character(treat_day_site_N2O$treatment)
treat_day_site_N2O$treatment <-
  as.factor(treat_day_site_N2O$treatment)
levels(treat_day_site_N2O$treatment)

treat_day_site_N2O$site <- as.character(treat_day_site_N2O$site)
treat_day_site_N2O$site <- as.factor(treat_day_site_N2O$site)
levels(treat_day_site_N2O$site)

#calculate standard errors
treat_day_site_N2O$SE <-
  treat_day_site_N2O$sd / sqrt(treat_day_site_N2O$n)
treat_day_site_N2O$SE[is.na(treat_day_site_N2O$SE)] <-
  0 #replacing NAs with 0
treat_day_site_N2O$SE

levels(treat_day_site_N2O$treatment)
levels(treat_day_site_N2O$site)

# convert date info in format 'mm/dd/yyyy'
treat_day_site_N2O$date <-
  paste(
    as.character(treat_day_site_N2O$month),
    as.character(treat_day_site_N2O$day),
    as.character(treat_day_site_N2O$year),
    sep = "/"
  )
treat_day_site_N2O$date
treat_day_site_N2O$date <-
  as.Date(treat_day_site_N2O$date, "%m/%d/%Y")

# convert date info in format 'mm/dd/yyyy'
treat_day_site_N2O_all$date <-
  paste(
    as.character(treat_day_site_N2O_all$month),
    as.character(treat_day_site_N2O_all$day),
    as.character(treat_day_site_N2O_all$year),
    sep = "/"
  )
treat_day_site_N2O_all$date
treat_day_site_N2O_all$date <-
  as.Date(treat_day_site_N2O_all$date, "%m/%d/%Y")

#formatting precipitation data
precip$Date <- as.Date(precip$Date, "%m/%d/%Y")

precip_summary <- precip %>%
  group_by(Date) %>%
  summarise(precip_mm = sum(ppt_cor)) %>%
  filter(Date >= "2017-02-08" & Date <= "2018-05-29")


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

#the within treat/day error is less than across measurement, so
#this is the plot we'll include with day/treatment combination summarized
#and faceted by treatment
cbPalette <-
  c("#009E73", "#56B4E9", "#D55E00", "#F0E442", "#0072B2") #"#CC79A7","#E69F00"

levels(treat_day_site_N2O$treatment)
levels(treat_day_site_N2O$treatment)[levels(treat_day_site_N2O$treatment) ==
                                       "80"] <- "maize"
levels(treat_day_site_N2O_all$treatment)
levels(treat_day_site_N2O_all$treatment)[levels(treat_day_site_N2O_all$treatment) ==
                                           "80"] <- "maize"

treat_day_site_N2O$treatment <-
  ordered(treat_day_site_N2O$treatment,
          levels = c("forest",  "soy", "maize"))

treat_day_site_N2O_all$treatment <-
  ordered(treat_day_site_N2O_all$treatment,
          levels = c("forest",  "soy", "maize"))


# facet label names
labs <- c("forest",  "soybean", "soybean-maize")
names(labs) <- c("forest",  "soy", "maize")

LU_N2O_plot <-
  ggplot() +
  geom_jitter(
    data = treat_day_site_N2O_all,
    aes(x = date, y = N2O_mg_N_m2_h),
    size = 3,
    width = 2,
    color = "darkgrey"
  ) +
  geom_point(
    data = treat_day_site_N2O,
    aes(x = date, y = mean_N2O_mg_N_m2_h,
        colour = treatment),
    shape = 17,
    size = 3
  ) +
  geom_line(data = treat_day_site_N2O, aes(
    x = date,
    y = mean_N2O_mg_N_m2_h,
    colour = factor(treatment)
  )) +
  geom_errorbar(
    data = treat_day_site_N2O,
    aes(
      x = date,
      ymax = mean_N2O_mg_N_m2_h + treat_day_site_N2O$SE,
      ymin = mean_N2O_mg_N_m2_h - treat_day_site_N2O$SE,
      colour = factor(treatment)
    ),
    na.rm = FALSE
  ) +
  ylab(expression(paste(
    N[2], "O flux (mg N ", " ", m ^ -2, " ",
    h ^ -1, ')', sep = ""
  ))) +
  xlab("Date") +
  theme_default() +
  geom_vline(data = treat_day_site_N2O,
             xintercept = min(treat_day_site_N2O$date) - 10) +
  geom_hline(data = treat_day_site_N2O, yintercept = 0) +
  ylim(-.25, 1.2) +
  scale_color_manual(
    name = "Treatment",
    values = c(cbPalette[1],  "blue", cbPalette[3]),
    breaks = c("forest", "soy", "maize"),
    labels = c("forest",  "soybean", "soybean-maize")
  ) +
  theme(legend.position = "bottom", legend.key = element_blank())  +
  facet_grid(~ treatment, labeller = labeller(treatment = labs)) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(breaks = as.Date(c(
    "2017-02-01", "2017-07-01",
    "2017-11-01",
    "2018-04-01"
  )),
  labels = c("Feb 17", "Jul 17", "Nov 17",
             "Apr 18"))

LU_N2O_plot


#linear mixed effects model on N2O fluxes ~ treatment and date
treat_day_site_N2O$treatment = factor(treat_day_site_N2O$treatment)
treat_day_site_N2O$date = factor(treat_day_site_N2O$date)
treat_day_site_N2O$site = factor(treat_day_site_N2O$site)

model = lmer(mean_N2O_mg_N_m2_h ~ treatment + date + (1 |
                                                        site), data = treat_day_site_N2O)
anova(model)  #date is significant, not treatment
lmerTest::rand(model)   #random effects  significant


#precipitation plot
precip_plot <- ggplot(data = precip_summary, aes(x = (Date),
                                                 y = precip_mm)) +
  geom_col(alpha = .8) +
  ylab("Precipitation (mm)") +
  xlab(" ") +
  theme_default() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_vline(xintercept = min(precip_summary$Date) - 10) +
  geom_hline(yintercept = 0)



# Calculating and plotting monthly fluxes -------------------------------------------------------------

#excluding harvest and planting treatments
dat <-
  filter(dat, LU == 1 &
           treatment != "planted" & treatment != "harvest")
dat <- filter(dat, site != "planted" & treatment != "harvest")

#removing those levels from treatment and site
dat$treatment <- as.character(dat$treatment)
dat$treatment <- as.factor(dat$treatment)
levels(dat$treatment)[levels(dat$treatment) == "80"] <- "maize"
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

#bootstrapping distributions of fluxes from each unique treatment, site, and sampling date to compile 5 flux observations for each sample event
set.seed(100)
sample_fluxes = NULL
flux_row_list <- list()
chamber_number <- 0
for (i in unique(chamb_lengths$unique_ID)) {
  chamber_number <- chamber_number + 1
  message(paste0("Chamber number:", chamber_number))
  tempdf = chamb_lengths[which(chamb_lengths$unique_ID == i),]
  site = rep(tempdf$site[1], 5)
  treatment = rep(tempdf$treatment[1], 5)
  julian_day = rep(tempdf$j_day_raw[1], 5)
  replicate = c("a", "b", "c", "d", "e")
  
  
  if (tempdf$n_chamb[1] == 5) {
    print(paste("At", i, ": 5 chambers"))
    flux_pool = sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 5, replace =
                           F)
  } else if (tempdf$n_chamb[1] == 4) {
    print(paste("At", i, ": 4 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 4, replace =
                             F),
                  sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 1, replace =
                             F))
  } else if (tempdf$n_chamb[1] == 3) {
    print(paste("At", i, ": 3 chambers"))
    flux_pool = c(sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 3, replace =
                             F),
                  sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 2, replace =
                             F))
  } else if (tempdf$n_chamb[1] == 2) {
    print(paste("At", i, ": 2 chambers"))
    flux_pool = c(
      sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 2, replace = F),
      sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 2, replace =
                 F),
      sample_n(as.data.frame(tempdf$N2O_mg_N_m2_h), 1, replace =
                 F)
    )
  }
  
  #creating randomly sampled pools of 5 flux measurements with some repetition for those with n<5
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
  tempdf <- tempdf[order(tempdf$julian_day),]
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

#calculate annual cumulative fluxes of N2O-N per hectare
cumulative_fluxes$kg_N2O_N_ha <-
  cumulative_fluxes$cumflux / 10 ^ 2 *
  365 / cumulative_fluxes$total_days #converts from mg N/m2/day to kg N/ha/year


# cumulative fluxes plots and stats ---------------------------------------

cumulative_N2O_fluxes <- cumulative_fluxes

#linear mixed effects model on cumulative fluxes ~ treatment
cumulative_N2O_fluxes$treatment = factor(cumulative_N2O_fluxes$treatment)

model2 = lmer(kg_N2O_N_ha ~ treatment + (1 |
                                           site), data = cumulative_N2O_fluxes)
anova(model2)  #date is significant, not treatment
emmeans(model2, specs = pairwise ~  "treatment")#post hoc test

#values for table s2
cumulative_fluxes %>%
  group_by(treatment) %>%
  summarise(mean = mean(kg_N2O_N_ha),
            SE = sd(kg_N2O_N_ha) / sqrt(15))

#levels order
cumulative_N2O_fluxes$site <- factor(
  cumulative_N2O_fluxes$site,
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

#reordering treatments
cumulative_N2O_fluxes$treatment <-
  ordered(cumulative_N2O_fluxes$treatment,
          levels = c("forest",  "soy", "maize"))

#plotting annual fluxes
annual_N2O_plot <-
  ggplot(data = cumulative_N2O_fluxes,
         aes(
           x = treatment,
           y = kg_N2O_N_ha,
           color = site,
           shape = site,
           group = treatment
         ),
  ) +
  geom_violin() +
  geom_jitter(size = 4, width = .05) +
  scale_shape_manual(name = "Site", values = (15:27)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.4) +
  #geom_segment(aes(x=0.1,xend=0.1,y=-5,yend=5), colour="black") +
  #geom_segment(aes(x=0,xend=3.5,y=0,yend=0),colour="black") +
  ylab(expression(
    paste(N[2], "O flux (kg ",  ~ N[2], "O-N ", " ", ha ^ -1, " ",
          yr ^ -1, ')', sep = "")
  )) +
  xlab(" ") +
  scale_color_manual(name = "Site", values = RColorBrewer::brewer.pal(8, 'Dark2')) +
  theme_default() +
  theme(legend.position = "bottom", legend.key = element_blank()) +
  scale_x_discrete(labels = c("forest",  "soybean", "soybean-maize"))

annual_N2O_plot
