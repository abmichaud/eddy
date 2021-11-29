# prework work
library(tidyverse)
library(lubridate)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(imputeTS)
library(dplyr)

# theme(axis.line = element_line(colour = "black"),
#panel.grid.major = element_blank(),
#panel.grid.minor = element_blank(),
#panel.border = element_blank(),
#panel.background = element_blank())
#is the same as theme_classic()

## Geochemical analysis ##

# Burrow wall image analysis #

burrow <- read_tsv("burrow.txt")

bur_ox <- ggplot(data=burrow) +
  geom_point(aes(x=day, y=area, color=mesocosm)) +
  scale_color_manual(name=NULL,
                     values=c("#B40F20", "#E58601", "#46ACC8"),
                     breaks=c("C", "D", "A"),
                     labels=c("48 uM", "106 uM", "280 uM")) +
  stat_smooth(method = lm, formula = y~x, aes(x=day, y=area, color=mesocosm), data = burrow) +
  stat_poly_eq(formula = y~x, data = burrow,
               aes(x=day, y=area, color=mesocosm, label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE, label.x = c(0.05, 0.05, 0.05), label.y = c(0.74, 0.98, 0.87)) +
  stat_cor(aes(x=day, y=area, color=mesocosm), label.x = c(0.5, 0.5, 0.5), label.y = c(4.8, 6.8, 5.9)) +
  labs(x="Time (d)",
       y="Fe-Oxide area surrounding burrow (cm2)") +
  theme_classic()

pdf("burrow_feox.pdf", height = 4, width = 6)
ggarrange(bur_ox, ncol = 1, nrow = 1)
dev.off()

# Porewater Fe #

# meso 1, with nereis, porewater Fe(II)
m1n.pwfe<-read_tsv("meso_one_nereis_pwfe_timecourse.txt")

estimate_Npwfe_by_date <- function(target_date, target_depth) {
  data_for_date <- m1n.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Npwfe_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Npwfe_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Npwfe_interp_depth<-crossing(tibble(date = unique(m1n.pwfe$date)), 
                           tibble(depth = seq(-6.5, 12.5, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Npwfe_by_date(date[1], depth))

Npwfe1_raster <- crossing(tibble(date = seq(ymd("2017-05-13"), ymd("2017-05-23"), by=0.2)),
                       tibble(depth = unique(Npwfe_interp_depth$depth))) %>%
  #filter(depth < -1) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Npwfe_by_depth(depth[1], date))

# meso 1, CONTROL, porewater Fe(II)
m1c.pwfe<-read_tsv("meso_one_cont_pwfe_timecourse.txt")

estimate_Cpwfe_by_date <- function(target_date, target_depth) {
  data_for_date <- m1c.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Cpwfe_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Cpwfe_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Cpwfe_interp_depth<-crossing(tibble(date = unique(m1c.pwfe$date)), 
                             tibble(depth = seq(-6.5, 11.5, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Cpwfe_by_date(date[1], depth))

Cpwfe1_raster <- crossing(tibble(date = seq(ymd("2017-05-13"), ymd("2017-05-23"), by=0.2)),
                         tibble(depth = unique(Cpwfe_interp_depth$depth))) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Cpwfe_by_depth(depth[1], date))

# meso 3, with nereis, porewater Fe(II)
m3n.pwfe<-read_tsv("meso_three_nereis_pwfe_timecourse.txt")

estimate_Npwfe3_by_date <- function(target_date, target_depth) {
  data_for_date <- m3n.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Npwfe3_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Npwfe3_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Npwfe3_interp_depth<-crossing(tibble(date = unique(m3n.pwfe$date)), 
                             tibble(depth = seq(-4, 14, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Npwfe3_by_date(date[1], depth))

Npwfe3_raster <- crossing(tibble(date = seq(ymd("2017-11-09"), ymd("2017-11-19"), by=0.2)),
                          tibble(depth = unique(Npwfe3_interp_depth$depth))) %>%
  #filter(depth < -1) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Npwfe3_by_depth(depth[1], date))

# meso 3, CONTROL, porewater Fe(II)
m3c.pwfe<-read_tsv("meso_three_cont_pwfe_timecourse.txt")

estimate_Cpwfe3_by_date <- function(target_date, target_depth) {
  data_for_date <- m3c.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Cpwfe3_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Cpwfe3_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Cpwfe3_interp_depth<-crossing(tibble(date = unique(m3c.pwfe$date)), 
                             tibble(depth = seq(-3, 15, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Cpwfe3_by_date(date[1], depth))

Cpwfe3_raster <- crossing(tibble(date = seq(ymd("2017-11-09"), ymd("2017-11-19"), by=0.2)),
                         tibble(depth = unique(Cpwfe3_interp_depth$depth))) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Cpwfe3_by_depth(depth[1], date))

# meso 4, with nereis, porewater Fe(II)
m4n.pwfe<-read_tsv("meso_four_nereis_pwfe_timecourse.txt")

estimate_Npwfe4_by_date <- function(target_date, target_depth) {
  data_for_date <- m4n.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Npwfe4_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Npwfe4_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Npwfe4_interp_depth<-crossing(tibble(date = unique(m4n.pwfe$date)), 
                             tibble(depth = seq(-5, 13, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Npwfe4_by_date(date[1], depth))

Npwfe4_raster <- crossing(tibble(date = seq(ymd("2018-06-10"), ymd("2018-06-20"), by=0.2)),
                          tibble(depth = unique(Npwfe4_interp_depth$depth))) %>%
  #filter(depth < -1) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Npwfe4_by_depth(depth[1], date))

# meso 4, CONTROL, porewater Fe(II)
m4c.pwfe<-read_tsv("meso_four_cont_pwfe_timecourse.txt")

estimate_Cpwfe4_by_date <- function(target_date, target_depth) {
  data_for_date <- m4c.pwfe %>%
    filter(date == target_date) %>%
    arrange(depth)
  approx(data_for_date$depth, data_for_date$pwfe, xout = target_depth)$y
}

estimate_Cpwfe4_by_depth <- function(target_depth, target_date) {
  data_for_depth <- Cpwfe4_interp_depth %>%
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$pwfe, xout = target_date)$y
}

Cpwfe4_interp_depth<-crossing(tibble(date = unique(m4c.pwfe$date)), 
                             tibble(depth = seq(-5, 13, length.out = 100))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Cpwfe4_by_date(date[1], depth))

Cpwfe4_raster <- crossing(tibble(date = seq(ymd("2018-06-10"), ymd("2018-06-20"), by=0.2)),
                          tibble(depth = unique(Cpwfe4_interp_depth$depth))) %>%
  group_by(depth) %>%
  mutate(pwfe = estimate_Cpwfe4_by_depth(depth[1], date))

# Plotting porewater Fe data
oxy_legend<-expression(paste("Oxygen (mL ",L^-1,")"))
fe_legend<-expression(paste("Fe(II) (",mu,"M)"))
Npwfe1_grid<-ggplot(Npwfe1_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m1n.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2017-05-12", "2017-05-24"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))
# another colorblind-friendly option: scale_fill_gradient2(low = "#009E73", mid = "#F0E442", high = "#D55E00", midpoint =

Cpwfe1_grid<-ggplot(Cpwfe1_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m1c.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2017-05-12", "2017-05-24"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

Npwfe3_grid<-ggplot(Npwfe3_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m3n.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2017-11-08", "2017-11-20"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

Cpwfe3_grid<-ggplot(Cpwfe3_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m3c.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2017-11-08", "2017-11-20"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

Npwfe4_grid<-ggplot(Npwfe4_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m4n.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2018-06-09", "2018-06-21"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

Cpwfe4_grid<-ggplot(Cpwfe4_raster, aes(date, depth, fill = pwfe)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "#F0E442", high = "#D55E00", midpoint = 220, name = fe_legend) +
  geom_point(data = m4c.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2018-06-09", "2018-06-21"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

pdf("meso_pwfe_time.pdf", height = 7, width = 10)
ggarrange(Npwfe1_grid, Cpwfe1_grid, Npwfe4_grid, Cpwfe4_grid, Npwfe3_grid, Cpwfe3_grid, ncol = 2, nrow = 3,
          labels = "AUTO", common.legend = TRUE, legend = "bottom")
dev.off()


## Fe(II) oxidation rate
Npwfe_interp_depth_equal<-crossing(tibble(date = unique(m1n.pwfe$date)), 
                             tibble(depth = seq(-6.5, 12.5, length.out = 190))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Npwfe_by_date(date[1], depth))

Cpwfe_interp_depth_equal<-crossing(tibble(date = unique(m1c.pwfe$date)), 
                             tibble(depth = seq(-6.5, 11.5, length.out = 180))) %>%
  group_by(date) %>%
  mutate(pwfe = estimate_Cpwfe_by_date(date[1], depth))

depth_equal <- round(as.numeric(Cpwfe_interp_depth_equal$depth), digits = 1)
Cpwfe1_depth_equal <- Cpwfe_interp_depth_equal %>% select(-depth)
Cpwfe1_depth_equal<-cbind(Cpwfe1_depth_equal, depth=depth_equal)

depth_equal <- round(as.numeric(Npwfe_interp_depth_equal$depth), digits = 1)
Npwfe1_depth_equal <- Npwfe_interp_depth_equal %>% select(-depth)
Npwfe1_depth_equal<-cbind(Npwfe1_depth_equal, depth=depth_equal) 

meso1 <- merge(Cpwfe1_depth_equal, Npwfe1_depth_equal, by = "depth", all = FALSE)
meso1_filt <- meso1 %>% filter(depth > 0) %>%
  mutate(feox = (pwfe.x - pwfe.y)) %>%
  filter(date.x == date.y)

ggplot(meso1_filt, aes(date.x, depth, fill = feox)) +
  geom_raster() +
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", name = fe_legend) +
  #geom_point(data = m1n.pwfe, aes(x = samp_date, y = samp_depth), size = 0.5) +
  scale_y_reverse() +
  ylim(15, -6) +
  xlab(NULL) +
  scale_x_date(limits = as.Date(c("2017-05-12", "2017-05-24"))) +
  ylab("Depth (cm)") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color = "white") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.position = c(0.947, 0.5), legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))



# Fe flux calculations
flux <- m1c.pwfe %>% filter(date == as.Date("2017-05-16"))
plot(flux$depth, flux$pwfe)
str(m1c.pwfe)

flux <- data.frame(mesocosm=numeric(), treatment=character(), flux=numeric())

x <- m1c.pwfe %>% filter(date == as.Date("2017-05-18"))

#z is data file
#mc is mesocosm ID
#tm is a treatment ID
#n is the days past start where the calculation should occur
#dif is the diffusion coefficient
#por is the porosity
#z is a tibble with date, depth, and porewater Fe2+ concentration in umoles per L
#on return: tibble with mesocosm ID, treatment ID, and Fe2+ flux umoles per cm-2 s-1

fitpwfe <- function(z, mc=1, tm, n=4, dif=0.00000582, por=0.85) {
  df <- data.frame(mesocosm=numeric(), treatment=numeric(), flux=numeric())
  startdate <- max(head(z$date))
  startdate <- startdate + n
  p1 <- filter(z, date == as.Date(startdate), depth > 0) 
  if(p1[1,3]>0) {
    slope <- z %>% filter(date == as.Date(startdate)) %>%
      filter(depth > 0, depth < 5) %>% lm(pwfe ~ depth, data = .)
    df %>% add_row(mesocosm = 1, 
                   flux = (slope[["coefficients"]][["depth"]]*dif*por))
  } else {
    df %>% add_row(mesocosm = 1, flux = 0)
  }
}

# @param x at least two columns, one depth and one measured value
# @param depth is numeric and deeper than 0 is positive
# @param var is the name of the measured value column
# @return a tibble with an added row of interpolated value at value depth

interp_at_depth <- function(x, depth=0, var="concentration") {
  p1 <- filter(x, depth > 0)
  if(nrow(p1)==0) return(x)
  p1 <- p1 %>% arrange(depth)
  if(p1$depth[1]>0) {
    slope <- x %>% filter(depth > 0, depth < 5) %>% lm(pwfe ~ depth, data = .)
    df %>% add_row(mesocosm = 1, 
                   flux = (slope[["coefficients"]][["depth"]]*dif*por))
  } else {
    df %>% add_row(mesocosm = 1, flux = 0)
  }
}

df <- data.frame(mesocosm=numeric(), flux=numeric())
startdate <- max(head(Npwfe1_raster$date))
startdate <- startdate + 4
p1 <- filter(Npwfe1_raster, date == as.Date(startdate), depth > 0) 
if(p1[1,3]>0) {
  slope <- Npwfe1_raster %>% filter(date == as.Date(startdate)) %>%
    filter(depth > 0, depth < 5) %>% lm(pwfe ~ depth, data = .)
  df %>% add_row(mesocosm = 1, 
                 flux = (slope[["coefficients"]][["depth"]]*0.00000582*0.85))
} else {
  df %>% add_row(mesocosm = 1, flux = 0)
  }

write.csv(fitpwfe(Cpwfe1_raster, 1, 0, 4, 0.00000582, 0.85), "m1c.csv")
write.csv(fitpwfe(Npwfe1_raster, 1, 25, 4, 0.00000582, 0.85), "m1n.csv")
m3c.pwfe <- m3c.pwfe %>% select(-c(timepoint, exp_type))
write.csv(fitpwfe(Cpwfe3_raster, 3, 0, 4, 0.00000582, 0.85), "m3c.csv")
m3n.pwfe <- m3n.pwfe %>% select(-c(timepoint, exp_type))
write.csv(fitpwfe(Npwfe3_raster, 3, 25, 4, 0.00000582, 0.85), "m3n.csv")
m4c.pwfe <- m4c.pwfe %>% select(-c(timepoint, exp_type))
write.csv(fitpwfe(Cpwfe4_raster, 4, 0, 4, 0.00000582, 0.85), "m4c.csv")
m4n.pwfe <- m4n.pwfe %>% select(-c(timepoint, exp_type))
write.csv(fitpwfe(Npwfe4_raster, 4, 25, 4, 0.00000582, 0.85), "m4n.csv")

# Poorly Crystalline particulate Fe and AVS content #
int_fes<-read_tsv("pcfe.txt")

pcfe <- ggplot(int_fes) +
  geom_boxplot(aes(x=meso, y=int_pcfe, color=treatment), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=int_pcfe, color=treatment)) +
  scale_color_manual(name=NULL,
                     values=c("black", "red", "blue"),
                     breaks=c("sieve", "nereis", "control"),
                     labels=c("Initial", "Bioturbated", "Control")) +
  scale_x_discrete(limits=c("initial", "one", "four", "three"),
                   labels=c("Initial Material", "280 umol", "106 umol", "48 umol")) +
  labs(x="Mesocosm", y="Poorly Crystalline Fe (umol Fe cm-2)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

pfes <- ggplot(int_fes) +
  geom_boxplot(aes(x=meso, y=int_avs, color=treatment), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=int_avs, color=treatment)) +
  scale_color_manual(name=NULL,
                     values=c("black", "red", "blue"),
                     breaks=c("sieve", "nereis", "control"),
                     labels=c("Initial", "Bioturbated", "Control")) +
  scale_x_discrete(limits=c("initial", "one", "four", "three"),
                   labels=c("Initial Material", "280 umol", "106 umol", "48 umol")) +
  labs(x="Mesocosm", y="Acid Volatile Sulfide (umol S cm-2)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

fe_flux <- read_tsv("fe_fluxes.txt")

flux_plot <- ggplot(fe_flux) +
  geom_bar(aes(x=as.factor(O2_conc), y=avg_flux, color=treatment, fill=treatment), position = "dodge", stat = "identity") +
  geom_errorbar(aes(x=as.factor(O2_conc), ymin=avg_flux-sd_flux, ymax=avg_flux+sd_flux, fill=treatment), 
                width=0.4, alpha=0.9, size=0.7, position=position_dodge(.9)) +
  scale_fill_manual(name=NULL,
                     values=c("red", "blue"),
                     breaks=c("Nereis", "control"),
                     labels=c("Bioturbated", "Control")) +
  scale_color_manual(name=NULL,
                    values=c("red", "blue"),
                    breaks=c("Nereis", "control"),
                    labels=c("Bioturbated", "Control")) +
  scale_x_discrete(limits=c("280", "106", "48"),
                   labels=c("280 umol", "106 umol", "48 umol")) +
  labs(x="Mesocosm", y="Fe(II) flux (umol cm-2 d-1)") +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"))

pdf("meso_solid_fes.pdf", height = 4, width = 14)
ggarrange(pcfe, pfes, flux_plot, ncol = 3, nrow = 1,
          labels = "AUTO")
dev.off()

# Mesocosm oxygen supplemental figure #
o2meso1 <- read_tsv("meso1_oxygen_measurements_17apr2020.txt")
o2meso3 <- read_tsv("meso3_oxygen_measurements_17apr2020.txt")
o2meso4 <- read_tsv("meso4_oxygen_measurements_17apr2020.txt")

o2meso1 <- mutate(o2meso1, time = row_number()) %>%
  mutate(day = time/1440) %>%
  select(-time)
o2meso3 <- mutate(o2meso3, time = row_number()) %>%
  mutate(day = time/86400) %>%
  select(-time)
o2meso4 <- mutate(o2meso4, time = row_number()) %>%
  mutate(day = time/86400) %>%
  select(-time)

oxy_plot <- ggplot() +
  geom_point(data=o2meso1, aes(x=day, y=treat), size = 1, shape = 17, color = "black") +
  geom_point(data=o2meso1, aes(x=day, y=cont), size = 1, shape = 18, color = "red") +
  geom_point(data=o2meso3, aes(x=day, y=treat), size = 1, shape = 17, color = "black") +
  geom_point(data=o2meso3, aes(x=day, y=cont), size = 1, shape = 18, color = "red") +
  geom_point(data=o2meso4, aes(x=day, y=treat), size = 1, shape = 17, color = "black") +
  geom_point(data=o2meso4, aes(x=day, y=cont), size = 1, shape = 18, color = "red") +
  labs(x="Time (d)", y="O2 (umol L-1)") +
  ylim(c(0,300)) +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), 
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0))

pdf("meso_oxy.pdf", height = 4, width = 6)
ggarrange(oxy_plot, ncol=1, nrow=1)
dev.off()

## Microbial Analysis ##
metadata<-read_tsv("eddy_meta.txt")
eddy.tax <- read_tsv("stability.final.pick.opti_mcc.0.15.cons.taxonomy") %>% 
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ";$", replacement = "")) %>%
  mutate(taxonomy=str_replace_all(string = taxonomy, pattern = ".__", replacement = "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

# @input mothur formatted and name shortened shared file (# of OTUS in each sample)
# @return tidy dataframe in long format with OTUs in % relative abundance
otu_table <- read_tsv("stability.final.pick.opti_mcc.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group)

rel_abun_rowSums <- function(x, y) {
  mutate(y, n=rowSums(x))
}

total_reads <- otu_table %>% 
  group_by(sample) %>% 
  group_map(rel_abun_rowSums) %>%
  bind_rows()

otu_table <- otu_table %>% mutate(reads=total_reads$n)

eddy.otu <- otu_table %>%
  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/reads)*100)

## Remove singletons and calculate relative abundances, use this below.
##@ return, get a long OTU table with relative abundances calculated with reads of each sample with singletons removed.
eddy.otu.single <- read_tsv("stability.final.pick.opti_mcc.shared", col_types = cols(Group=col_character())) %>%
  select(-label, -numOtus)
sample <- eddy.otu.single$Group
eddy.otu.single <- select(eddy.otu.single, -Group)
eddy.otu.single <- eddy.otu.single[,colSums(eddy.otu.single)>1]
eddy.otu.single <- cbind(eddy.otu.single, sample)

rel_abun_rowSums <- function(x, y) {
  mutate(y, n=rowSums(x))
}

total_reads_single <- eddy.otu.single %>% 
  group_by(sample) %>% 
  group_map(rel_abun_rowSums) %>%
  bind_rows()

eddy.otu.single <- eddy.otu.single %>% mutate(reads=total_reads_single$n)

eddy.otu.single <- eddy.otu.single %>%
  pivot_longer(cols = c(-sample, -reads), names_to = "otu", values_to = "count") %>%
  mutate(rel_abund=(count/reads)*100)


# @input reformatted OTU table and OTU taxonomy table
# @return % rel abund and taxonomy for every OTU from every sample
otu.tax <- inner_join(eddy.otu.single, eddy.tax)

# @input metadata for all samples and % rel abund and taxonomy for every OTU from every sample
# @return complete dataframe of all OTUs from every sample with everything we want to know about each OTU
otu.tax.meta <- inner_join(otu.tax, metadata, by=c("sample"="sample"))

gal<-filter(otu.tax.meta, family%in%c("Gallionellaceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "sample"))

gal_plot<-ggplot(gal) +
  geom_boxplot(aes(x=meso, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=rel_abund, color=treat)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("initial", "nereis", "cont"),
                     labels=c("Initial Mud", "+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("eddy", "one", "four", "three"),
                   labels=c("Eddy Mud", "280 ",mu,"M", "106 uM", "48 uM")) +
  labs(title="Gallionellaceae", x="Mesocosm",
       y="Relative Abundance (%)") +
  theme_classic()




feox<-filter(otu.tax.meta, family%in%c("Mariprofundaceae", "Gallionellaceae")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "sample"))

feox_plot<-ggplot(feox) +
  geom_boxplot(aes(x=meso, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=rel_abund, color=treat)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("initial", "nereis", "cont"),
                     labels=c("Initial Mud", "+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("eddy", "one", "three", "four"),
                   labels=c("Eddy Mud", "Meso1", "Meso3", "Meso4")) +
  labs(title="Zetaproteobacteria", x="Mesocosm",
       y="Relative Abundance (%)") +
  theme_classic()

feox_plot <- ggplot(feox) +
  geom_boxplot(aes(x=meso, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=rel_abund, color=treat)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("initial", "nereis", "cont"),
                     labels=c("Initial Mud", "+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("eddy", "one", "four", "three"),
                   labels=c("Eddy Mud", "280 uM", "106 uM", "48 uM")) +
  scale_shape_manual(values=c(17, 18, 19),
                     breaks = c("initial", "nereis", "cont"),
                     labels=c("initial", "nereis", "control")) +
  labs(y="Relative Abundance (%)", x=NULL) +
  ylim(0.2,0.8) +
  theme_classic()

cable<-filter(otu.tax.meta, family%in%"Desulfobulbaceae") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "sample"))

cable_plot<-ggplot(cable) +
  geom_boxplot(aes(x=meso, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=rel_abund, color=treat)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("initial", "nereis", "cont"),
                     labels=c("Initial Mud", "+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("eddy", "one", "four", "three"),
                   labels=c("Eddy Mud", "280 uM", "106 uM", "48 uM")) +
  labs(x="Mesocosm Oxygen Concentration", y=NULL) +
  theme_classic()

desulfo <-filter(otu.tax.meta, family%in%"Desulfobacteraceae") %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "sample"))

desulfo_plot <- ggplot(desulfo) +
  geom_boxplot(aes(x=meso, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(shape=19, size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=meso, y=rel_abund, color=treat)) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("initial", "nereis", "cont"),
                     labels=c("Initial Mud", "+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("eddy", "one", "four", "three"),
                   labels=c("Eddy Mud", "280 uM", "106 uM", "48 uM")) +
  labs(x=NULL, y=NULL) +
  theme_classic()

pdf("eddy_taxa.pdf", height = 4, width = 12)
ggarrange(feox_plot, cable_plot, desulfo_plot, ncol=3, nrow=1, 
          labels = "AUTO", common.legend = TRUE, legend = "right")
dev.off()


sulfox <-filter(otu.tax.meta, genus%in%c("Sulfurimonas", "Thiobacillus")) %>%
  group_by(sample) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  inner_join(., metadata, by=c("sample" = "sample")) %>%
  filter(meso != "eddy")

ggplot(sulfox) +
  geom_point(aes(x=depth, y=rel_abund, color=treat), outlier.shape = NA) +
  geom_jitter(size=2, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), 
              aes(x=depth, y=rel_abund, color=treat)) +
  #scale_color_manual(name=NULL,
                    # values=c("red", "black"),
                   #  breaks=c("nereis", "cont"),
                    # labels=c("+Nereis", "-Nereis")) +
  scale_x_discrete(limits=c("1", "3", "5", "10"),
                   labels=c("1 cm", "3 cm", "5 cm", "10 cm")) +
  ylim(0,0.4) +
  labs(x="Depth",
       y="Relative Abundance (%)") +
  theme_classic()

ggplot(sulfox) +
  geom_point(aes(x=rel_abund, y=depth, color=meso, shape=treat), outlier.shape = NA) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black", "green"),
                     c("eddy", "one", "four", "three"),
                     labels=c("Eddy Mud", "280uM O2", "48 uM O2", "106 uM O2")) +
  xlim(0,0.3) +
  scale_y_reverse() +
  labs(x="Relative Abundance (%)",
       y="Depth (cm)") +
  theme_classic()




