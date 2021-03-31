# Packages

library(tidyverse)
library(ggridges)
library(cowplot)
library(drc)
library(effsize)
library(reshape2)
library(grid)
library(car)
library(readxl)

# Files

subsets_BPS <- read_csv("~/Desktop/Research/Hot plate assay paper/Data analysis/BPS_subsets_filtered.csv")
RHP_BPS_control_raw_data <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/BPS control raw data files/RHP BPS control raw data.xlsx")
SHP_BPS_control_raw_data_ <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/BPS control raw data files/SHP BPS control raw data .xlsx")
BPS_SHP_THC_DR_data <- read_excel("~/Desktop//Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/Excel files/BPS SHP THC DR data.xlsx")
BPS_RHP_THC_DR <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/Excel files/BPS RHP THC DR.xlsx")
BPS_morphine_DR_data <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/Excel files/BPS morphine DR data.xlsx")
BPS_G_power_plots_copy <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/Excel files/BPS G_power plots copy.xlsx")

## Processing files

### Control subpopulations

subsets_BPS$Subset <- factor(subsets_BPS$Subset, levels = c("Untreated", "Vehicle-treated", "Combined"))
subsets_BPS$Protocol <- factor(subsets_BPS$Protocol, levels = c("Standard hot plate", "Ramped hot plate"))

subsets_BPS %>%
  group_by(Protocol, Subset) %>%
  summarize(count = n())

### SHP QQ-density plot

SHP_untreated_vehicle <- SHP_BPS_control_raw_data_ %>%
  gather(key = "Subset", value = "Latency", -Experiment, -rep_measurement) %>%
  filter(!is.na(Latency)) %>%
  mutate(latency_correction = (Latency - 9) / 111,
         log_Latency = log10(Latency),
         log_latency_correction = (log_Latency - log10(9)) / 1.125)

bartlett.test(Latency ~ Subset, data = SHP_untreated_vehicle)
t.test(Latency ~ Subset, data = SHP_untreated_vehicle, var.equal = F)
cohen.d(d = SHP_untreated_vehicle$Latency, f = SHP_untreated_vehicle$Subset, 
        hedges.correction = T, conf.level = 0.95)

SHP_combined_BPS <- SHP_untreated_vehicle %>%
  mutate(latency_correction = (Latency - 9) / 111,
         log_Latency = log10(Latency),
         log_latency_correction = (log_Latency - log10(9)) / 1.125)

SHP_normal_combined_slope <- 
  diff(quantile(SHP_combined_BPS$latency_correction, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

SHP_normal_combined_intercept <- 
  quantile(SHP_combined_BPS$latency_correction, 0.25, na.rm = T) - 
  SHP_normal_combined_slope * qnorm(0.25)

### RHP QQ-density plot

RHP_untreated_vehicle <- RHP_BPS_control_raw_data %>%
  gather(key = "Subset", value = "Latency", -Experiment, -rep_measurement) %>%
  filter(!is.na(Latency)) %>%
  mutate(latency_correction = (Latency - 131.19) / 94.08,
         log_Latency = log10(Latency),
         log_latency_correction = (log_Latency - log10(131.19)) / 0.235)

bartlett.test(Latency ~ Subset, data = RHP_untreated_vehicle)
t.test(Latency ~ Subset, data = RHP_untreated_vehicle, var.equal = T)
cohen.d(d = RHP_untreated_vehicle$Latency, f = RHP_untreated_vehicle$Subset, 
        hedges.correction = T, conf.level = 0.95)

RHP_combined_BPS <- RHP_untreated_vehicle %>%
  mutate(latency_correction = (Latency - 131.19) / 94.08,
         log_Latency = log10(Latency),
         log_latency_correction = (log_Latency - log10(131.19)) / 0.235)

### Dose-response studies: Reference compounds

#### Standard hot plate

##### Pure THC

logSHP_DR_scaled <- BPS_SHP_THC_DR_data %>% 
  dplyr::select(-Assay) %>%
  mutate(Effect_correction = (Effect - 9) / 111,
         log_Effect_correction = (log_Effect - 0.9542425) / 1.125)

leveneTest(log_Effect_correction ~ factor(Dose), data = logSHP_DR_scaled)
bartlett.test(log_Effect_correction ~ factor(Dose), data = logSHP_DR_scaled)

logSHP_DR_scaled %>% 
  group_by(Dose) %>%
  summarize(n = n())

##### Morphine

logSHP_morphine_DR <- BPS_morphine_DR_data %>%
  filter(Assay == "logSHP") %>%
  melt(id = c("Assay")) %>%
  filter(!is.na(value)) %>%
  mutate(log_latency_correction = (value - 0.9542425) / 1.125)

logSHP_morphine_DR$variable <- as.numeric(as.character(logSHP_morphine_DR$variable))

logSHP_morphine_DR <- logSHP_morphine_DR %>% dplyr::rename(Dose = variable)
logSHP_morphine_DR <- logSHP_morphine_DR %>% dplyr::rename(Latency = value)

leveneTest(log_latency_correction ~ factor(Dose), data = logSHP_morphine_DR)
bartlett.test(log_latency_correction ~ factor(Dose), data = logSHP_morphine_DR)

#### Ramped hot plate

##### Pure THC

logRHP_DR_scaled <- BPS_RHP_THC_DR %>%
  mutate(Effect_correction = (Effect - 131.19) / 93.790,
         log_Effect_correction = (log_Effect - 2.117901) / 0.235)

leveneTest(log_Effect_correction ~ factor(Dose), data = logRHP_DR_scaled)
bartlett.test(log_Effect_correction ~ factor(Dose), data = logRHP_DR_scaled)

##### Morphine

RHP_morphine_DR <- BPS_morphine_DR_data %>%
  filter(Assay == "logRHP") %>%
  melt(id = c("Assay")) %>%
  filter(!is.na(value)) %>%
  mutate(log_latency_correction = (value - 2.117901) / 0.235)

leveneTest(log_latency_correction ~ factor(Dose), data = RHP_morphine_DR)
bartlett.test(log_latency_correction ~ factor(Dose), data = RHP_morphine_DR)

RHP_morphine_DR$variable <- as.numeric(as.character(RHP_morphine_DR$variable))

RHP_morphine_DR <- RHP_morphine_DR %>% dplyr::rename(Dose = variable)
RHP_morphine_DR <- RHP_morphine_DR %>% dplyr::rename(Latency = value)

### ANOVA

#### SHP

##### THC

logSHP_DR_scaled

logSHP_THC_DR_spread <- logSHP_DR_scaled %>%
  dplyr::select(Dose, log_Effect_correction) %>%
  dplyr::mutate(Dose = factor(Dose, ordered = T),
                row_id=1:n()) %>%
  tidyr::spread(key = Dose, value = log_Effect_correction) %>%
  dplyr::select(-row_id) 

write.table(logSHP_THC_DR_spread, "./logSHP_THC_DR_spread.txt", sep="\t", row.names = F)

##### Morphine

logSHP_morphine_DR

logSHP_morphine_DR_spread <- logSHP_morphine_DR %>%
  dplyr::select(Dose, log_latency_correction) %>%
  dplyr::mutate(Dose = factor(Dose, ordered = T),
                row_id=1:n()) %>%
  tidyr::spread(key = Dose, value = log_latency_correction) %>%
  dplyr::select(-row_id) 

write.table(logSHP_morphine_DR_spread, "./logSHP_morphine_DR_spread.txt", sep="\t", row.names = F)

#### RHP

##### THC

logRHP_DR_scaled

logRHP_THC_DR_spread <- logRHP_DR_scaled %>%
  dplyr::select(Dose, log_Effect_correction) %>%
  dplyr::mutate(Dose = factor(Dose, ordered = T),
                row_id=1:n()) %>%
  tidyr::spread(key = Dose, value = log_Effect_correction) %>%
  dplyr::select(-row_id) 

write.table(logRHP_THC_DR_spread, "./logRHP_THC_DR_spread.txt", sep="\t", row.names = F)

##### Morphine

RHP_morphine_DR

logRHP_morphine_DR_spread <- RHP_morphine_DR %>%
  dplyr::select(Dose, log_latency_correction) %>%
  dplyr::mutate(Dose = factor(Dose, ordered = T),
                row_id=1:n()) %>%
  tidyr::spread(key = Dose, value = log_latency_correction) %>%
  dplyr::select(-row_id) 

write.table(logRHP_morphine_DR_spread, "./logRHP_morphine_DR_spread.txt", sep="\t", row.names = F)


### G*Power

#### Alpha = 0.05, 0.01 SHP

G_power <- BPS_G_power_plots_copy
G_power$alpha <- factor(G_power$alpha, ordered = T, levels = c(0.05, 0.01))
G_power$Treatment <- factor(G_power$Treatment, ordered = T, levels = c("Morphine", "THC"))

#### Combined G*Power plot

x_1 <- rep(1, 455)
x_2 <- rep(2, 455)
x_3 <- rep(3, 455)
x_4 <- rep(4, 455)

binded_x_ <- cbind(x_1,
                   x_2,
                   x_3,
                   x_4) 

binded_df <- as.data.frame(binded_x_)
melted_df <- melt(binded_df) %>% 
  dplyr::select(-variable)
melted_df$value <- factor(melted_df$value, levels = c(1, 2, 3, 4))

Combined_power_plot <- G_power %>%
  filter(alpha == 0.05) %>%
  dplyr::select(-SHP, -RHP, -alpha) %>%
  gather(key = "Protocol", value = "Value", -Sample_size, -Groups, -Treatment) %>%
  separate(Protocol, into = c("empty", "Protocol"), sep = "log") %>%
  dplyr::select(-empty) %>%
  mutate(Protocol = gsub("SHP", "Standard\nhot plate", Protocol),
         Protocol = gsub("RHP", "Ramped\nhot plate", Protocol))

Combined_power_plot$Protocol <- factor(Combined_power_plot$Protocol,
                                       levels = c("Standard\nhot plate", "Ramped\nhot plate"))

combined_melt <- cbind(Combined_power_plot, melted_df)
Power_plot <- as.data.frame(combined_melt)
Power_plot_2 <- Power_plot %>%
  mutate(line_SHP_THC = ifelse(value == 1, 32, NA),
         line_RHP_THC = ifelse(value == 2, 28, NA),
         line_SHP_morphine = ifelse(value == 3, 23, NA),
         line_RHP_morphine = ifelse(value == 4, 25, NA)) %>%
  gather(key = "factored_lines", value = "FLine",
         -Sample_size, -Groups, -Treatment, -Protocol, -Value, -value)

Power_plot_2$FLine <- factor(Power_plot_2$FLine, levels = c("32", "28", "23", "25"))

vline_df <- data.frame(z = levels(Power_plot_2$Treatment),
                       vl = c(32, 28, 23, 25))

##### Grouping data in order to allow for multiple lines specific to each panel

mean.data <- aggregate(x = Power_plot_2$Value, # use the y values
                       by = Power_plot_2[c("Protocol", "Treatment")], # group by Group and then by Subgroup
                       FUN = function(x) {
                         signif(mean(x), 4) # calculate the mean keep 4 significant numbers
                       }
)

colnames(mean.data) <- c("Protocol", "Treatment", "Average")

# Figures

## Control subpopulations

ggplot(subsets_BPS, aes(x = Subset, y = Latency)) +
  geom_point(aes(colour = Protocol), position = position_jitter(width = 0.20), 
             alpha = 0.25, size = 3, pch = 19) +
  facet_wrap(~ Protocol, scales = "free_y") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", size = 1, width = 0.5, alpha = 0.75) +
  stat_summary(fun.y = mean, geom = "point", pch = 15, size = 4, alpha = 0.85) +
  scale_colour_manual(values = c("#000000", "#404040"), guide = F) +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.x = element_text(size = 14, colour = "black", family = "Arial"),
        legend.text = element_text(size = 14, family = "Arial"),
        panel.background = element_rect(),
        panel.grid.major = element_line(colour = "grey80", size = 0.3),
        panel.grid.minor = element_line(colour = "grey88", size = 0.25),
        strip.background = element_blank()) +
  labs(x = " ",y = "Latency (s)")

### Image: Control subpopulation
ggsave(filename = "control_subpopulations.tiff", width = 10, height = 5)

## SHP QQ-density plot

### Untransformed

SHP_normal_combined_panel1 <- ggplot(SHP_combined_BPS, aes(sample = latency_correction)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = SHP_normal_combined_slope, intercept = SHP_normal_combined_intercept),
              size = 1, col = "black", alpha = 0.5, lty = 1) +
  labs(y = "Proportion of effect", x = "Theoretical distribution",
       title = "Standard hot plate (untreated + vehicle)",
       subtitle = "Quantile-Quantile and Density plots")  +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text = element_text(family = "Arial", size = 14)) +
  annotate(geom = "text", label = "n = 108", x = -2, y = 0.75, 
           size = 5, family = "Arial")

SHP_normal_combined_dens <- axis_canvas(SHP_normal_combined_panel1, axis = "y") +
  geom_vridgeline(data = SHP_combined_BPS, aes(y = latency_correction, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.5, size = 0.75, trim = F, fill = "#000000") 

SHP_normal_combined_pp1 <- insert_yaxis_grob(SHP_normal_combined_panel1, SHP_normal_combined_dens, grid::unit(0.2, "null"),
                                             position = "right")

ggdraw(SHP_normal_combined_pp1)

### Log-transformed

SHP_log_combined_slope <- 
  diff(quantile(SHP_combined_BPS$log_latency_correction, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

SHP_log_combined_intercept <- 
  quantile(SHP_combined_BPS$log_latency_correction, 0.25, na.rm = T) - 
  SHP_log_combined_slope * qnorm(0.25)

SHP_log_combined_panel1 <- ggplot(SHP_combined_BPS, aes(sample = log_latency_correction)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = SHP_log_combined_slope, intercept = SHP_log_combined_intercept),
              size = 1, col = "black", alpha = 0.5, lty = 1) +
  labs(y = "Proportion of log-effect", x = "Theoretical distribution",
       title = "Standard hot plate (untreated + vehicle)",
       subtitle = "Log-transform: Quantile-Quantile and Density plots")  +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, .2)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text = element_text(family = "Arial", size = 14)) +
  annotate(geom = "text", label = "n = 108", x = -2, y = 0.75, 
           size = 5, family = "Arial")

SHP_log_combined_dens <- axis_canvas(SHP_log_combined_panel1, axis = "y") +
  geom_vridgeline(data = SHP_combined_BPS, aes(y = log_latency_correction, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.5, size = 0.75, trim = F, fill = "#000000") 

SHP_log_combined_pp1 <- insert_yaxis_grob(SHP_log_combined_panel1, SHP_log_combined_dens, grid::unit(0.2, "null"),
                                          position = "right")

ggdraw(SHP_log_combined_pp1)

### Image: SHP QQ-density plot
SHP_QQ_density_plot <- plot_grid(SHP_normal_combined_pp1, SHP_log_combined_pp1, align = "h")
ggsave(filename = "SHP_QQ_density_plot.tiff", plot = SHP_QQ_density_plot, width = 12, height = 5)

## RHP QQ-density plot

### Untransformed 

RHP_normal_combined_slope <- 
  diff(quantile(RHP_combined_BPS$latency_correction, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

RHP_normal_combined_intercept <- 
  quantile(RHP_combined_BPS$latency_correction, 0.25, na.rm = T) - 
  RHP_normal_combined_slope * qnorm(0.25)

RHP_normal_combined_panel1 <- ggplot(RHP_combined_BPS, aes(sample = latency_correction)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = RHP_normal_combined_slope, intercept = RHP_normal_combined_intercept),
              size = 1, col = "black", alpha = 0.5, lty = 1) +
  labs(y = "Proportion of effect", x = "Theoretical distribution",
       title = "Ramped hot plate (untreated + vehicle)",
       subtitle = "Quantile-Quantile and Density plots")  +
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, .1)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text = element_text(family = "Arial", size = 14)) +
  annotate(geom = "text", label = "n = 910", x = -2, y = 0.7, 
           size = 5, family = "Arial")

RHP_normal_combined_dens <- axis_canvas(RHP_normal_combined_panel1, axis = "y") +
  geom_vridgeline(data = RHP_combined_BPS, aes(y = latency_correction, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.5, size = 0.75, trim = F, fill = "#404040") 

RHP_normal_combined_pp1 <- insert_yaxis_grob(RHP_normal_combined_panel1, RHP_normal_combined_dens, grid::unit(0.2, "null"),
                                             position = "right")

ggdraw(RHP_normal_combined_pp1)

### Log-transformed

RHP_log_combined_slope <- 
  diff(quantile(RHP_combined_BPS$log_latency_correction, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

RHP_log_combined_intercept <- 
  quantile(RHP_combined_BPS$log_latency_correction, 0.25, na.rm = T) - 
  RHP_log_combined_slope * qnorm(0.25)

RHP_log_combined_panel1 <- ggplot(RHP_combined_BPS, aes(sample = log_latency_correction)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = RHP_log_combined_slope, intercept = RHP_log_combined_intercept),
              size = 1, col = "black", alpha = 0.5, lty = 1) +
  labs(y = "Proportion of log-effect", x = "Theoretical distribution",
       title = "Ramped hot plate (untreated + vehicle)",
       subtitle = "Log-transform: Quantile-Quantile and Density plots")  +
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, .1)) +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 14), 
        axis.text = element_text(family = "Arial", size = 14)) +
  annotate(geom = "text", label = "n = 910", x = -2, y = 0.7, 
           size = 5, family = "Arial")

RHP_log_combined_dens <- axis_canvas(RHP_log_combined_panel1, axis = "y") +
  geom_vridgeline(data = RHP_combined_BPS, aes(y = log_latency_correction, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.5, size = 0.75, trim = F, fill = "#404040") 

RHP_log_combined_pp1 <- insert_yaxis_grob(RHP_log_combined_panel1, RHP_log_combined_dens, grid::unit(0.2, "null"),
                                          position = "right")

ggdraw(RHP_log_combined_pp1)

### Image: RHP QQ-density plot

RHP_QQ_density_plot <- plot_grid(RHP_normal_combined_pp1, RHP_log_combined_pp1, align = "h")
ggsave(filename = "RHP_QQ_density_plot.tiff", plot = RHP_QQ_density_plot, width = 12, height = 5)

## Dose response curves: reference compounds

### Standard hot plate

#### Pure THC

logSHP_scaled_fit <- drm(log_Effect_correction ~ Dose, data = logSHP_DR_scaled, 
                         fct = LL.4(fixed = c(NA, NA, 1, NA), 
                                    names = c("Slope","Lower Limit","Upper Limit","ED50")))
logSHP_scaled_line <- expand.grid(Dose = exp(seq(log(max(logSHP_DR_scaled$Dose)),
                                                 log(min(logSHP_DR_scaled$Dose)),length=100))) 

logSHP_scaled <- predict(logSHP_scaled_fit,newdata=logSHP_scaled_line,interval="confidence") 
logSHP_scaled_line$p <- logSHP_scaled[,1]
logSHP_scaled_line$pmin <- logSHP_scaled[,2]
logSHP_scaled_line$pmax <- logSHP_scaled[,3]

aov_logSHP_scaled <- logSHP_DR_scaled %>%
  spread(key = Dose, value = log_Effect_correction)

logTHC_SHP_scaled_graph <- ggplot(logSHP_DR_scaled, aes(x = Dose, y = log_Effect_correction)) +
  geom_point(colour = "black", fill = "black", alpha = 0.25, size = 2) + 
  geom_line(data = logSHP_scaled_line, aes(x = Dose,y = p)) + 
  theme_bw() +
  labs(title = "Standard hot plate: THC", 
       subtitle = "n = 4-20 per group", x = "Dose (mg/kg)", y = "Proportion of log-effect") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", colour = "black", width = 0.25) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "#000000", alpha = 0.75, width = 0.2, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", colour = "#000000", alpha = 0.85, size = 3, pch = 15) +
  geom_ribbon(data = logSHP_scaled_line, aes(x = Dose,y = p, ymin = pmin, ymax = pmax), alpha = 0.15) +
  scale_x_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10), 
                     labels =c("Vehicle", "-1", "0", "1")) +
  theme(text = element_text(family = "Arial", size = 14),
        axis.text = element_text(family = "Arial", size = 14)) +
  geom_abline(slope = 0, intercept = 1, lty = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  annotation_logticks(sides = "b")

ED(logSHP_scaled_fit, c(25, 50, 75), interval = "delta")

#### Morphine

logmorphine_SHPfit <- drm(log_latency_correction ~ Dose, data = logSHP_morphine_DR, 
                          fct = LL.4(fixed = c(NA, NA, 1, NA), 
                                     names = c("Slope","Lower Limit","Upper Limit","ED50")))
logmorphine_SHPline <- expand.grid(Dose = exp(seq(log(max(logSHP_morphine_DR$Dose)),
                                                  log(min(logSHP_morphine_DR$Dose)),length=100))) 
logmorphine_SHP <- predict(logmorphine_SHPfit,newdata=logmorphine_SHPline,interval="confidence") 
logmorphine_SHPline$p <- logmorphine_SHP[,1]
logmorphine_SHPline$pmin <- logmorphine_SHP[,2]
logmorphine_SHPline$pmax <- logmorphine_SHP[,3]

logmorphine_SHP_graph <- ggplot(logSHP_morphine_DR, aes(x = Dose, y = log_latency_correction)) +
  geom_point(colour = "black", fill = "black", alpha = 0.25, size = 3) + 
  geom_line(data = logmorphine_SHPline, aes(x = Dose,y = p)) + 
  theme_bw() +
  labs(title = "Standard hot plate: Morphine", 
       subtitle = "n = 7 per group", x = "Dose (mg/kg)", y = "Proportion of log-effect") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", colour = "black", width = 0.25) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "#000000", alpha = 0.75, width = 0.2, size = .75) +
  stat_summary(fun.y = mean, geom = "point", colour = "#000000", alpha = 0.85, size = 3, pch = 15) +
  geom_ribbon(data = logmorphine_SHPline, aes(x = Dose,y = p, ymin = pmin, ymax = pmax), alpha = 0.2) +
  scale_x_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10), 
                     labels =c("Vehicle", "-1", "0", "1")) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, family = "Arial")) +
  geom_abline(slope = 0, intercept = 1, lty = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1.24), breaks = seq(0, 1, 0.2)) +
  annotation_logticks(sides = "b")

ED(logmorphine_SHPfit, c(25, 50, 75), interval = "delta")

### Ramped hot plate

#### Pure THC

logRHPfit_scaled <- drm(log_Effect_correction ~ Dose, data = logRHP_DR_scaled, 
                        fct = LL.4(fixed = c(NA, NA, 1, NA), names = c("Slope","Lower Limit","Upper Limit","ED50")))
logRHPline_scaled <- expand.grid(Dose = exp(seq(log(max(logRHP_DR_scaled$Dose)),
                                                log(min(logRHP_DR_scaled$Dose)),length=100))) 
logRHP_scaled <- predict(logRHPfit_scaled,newdata=logRHPline_scaled,interval="confidence") 
logRHPline_scaled$p <- logRHP_scaled[,1]
logRHPline_scaled$pmin <- logRHP_scaled[,2]
logRHPline_scaled$pmax <- logRHP_scaled[,3]

logTHC_scaled_RHP_graph <- ggplot(logRHP_DR_scaled, aes(x = Dose, y = log_Effect_correction)) +
  geom_point(colour = "black", fill = "black", alpha = 0.25, size = 2) + 
  geom_line(data = logRHPline_scaled, aes(x = Dose,y = p)) + 
  theme_bw() +
  labs(title = "Ramped hot plate: THC",
       subtitle = "n = 8-16 per group", x = "Dose (mg/kg)", y = "Proportion of log-effect") +
  scale_x_continuous(trans = "log10", breaks = c(0.01, 0.1, 1.0, 10), 
                     labels =c("Vehicle", "-1", "0", "1")) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", colour = "black", width = 0.25) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "#404040", alpha = 0.75, width = 0.2, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", colour = "#404040", alpha = 0.85, size = 3, shape = 15) +
  geom_ribbon(data = logRHPline_scaled, aes(x = Dose,y = p, ymin = pmin, ymax = pmax), alpha = 0.15) +
  theme(text = element_text(family = "Arial", size = 14),
        axis.text = element_text(family = "Arial", size = 14),
        title = element_text(family = "Arial", size = 14)) +
  geom_abline(slope = 0, intercept = 1, lty = 3, alpha = 0.8, colour = "black") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  annotation_logticks(sides = "b")

ED(logRHPfit_scaled, c(25, 50, 75), interval = "delta")

#### Moprhine

morphine_RHPfit <- drm(log_latency_correction ~ Dose, data = RHP_morphine_DR, 
                       fct = LL.4(fixed = c(NA, NA, 1, NA), 
                                  names = c("Slope","Lower Limit","Upper Limit","ED50")))
morphine_RHPline <- expand.grid(Dose = exp(seq(log(max(RHP_morphine_DR$Dose)),
                                               log(min(RHP_morphine_DR$Dose)),length=100))) 
morphine_RHP <- predict(morphine_RHPfit,newdata=morphine_RHPline,interval="confidence") 
morphine_RHPline$p <- morphine_RHP[,1]
morphine_RHPline$pmin <- morphine_RHP[,2]
morphine_RHPline$pmax <- morphine_RHP[,3]

morphine_RHP_graph <- ggplot(RHP_morphine_DR, aes(x = Dose, y = log_latency_correction)) +
  geom_point(colour = "black", fill = "black", alpha = 0.25, size = 3) + 
  geom_line(data = morphine_RHPline, aes(x = Dose,y = p)) + 
  theme_bw() +
  labs(title = "Ramped hot plate: Morphine", 
       subtitle = "n = 8 per group", x = "Dose (mg/kg)", y = "Proportion of log-effect") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", colour = "black", width = 0.25) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "#404040", alpha = 0.75, width = 0.2, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", colour = "#404040", alpha = 0.85, size = 3, pch = 15) +
  geom_ribbon(data = morphine_RHPline, aes(x = Dose,y = p, ymin = pmin, ymax = pmax), alpha = 0.2) +
  scale_x_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10), 
                     labels =c("Vehicle", "-1", "0", "1")) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, family = "Arial")) +
  geom_abline(slope = 0, intercept = 1, lty = 3, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1.24), breaks = seq(0, 1, 0.2)) +
  annotation_logticks(sides = "b")

ED(morphine_RHPfit, c(25, 50, 75), interval = "delta")

### Images

#### Pure THC

THC_DR_plots <- plot_grid(logTHC_SHP_scaled_graph, logTHC_scaled_RHP_graph, align = "h")
ggsave(filename = "THC_DR_plots.tiff", plot = THC_DR_plots, width = 10, height = 5)

#### Morphine

Morphine_DR_plots <- plot_grid(logmorphine_SHP_graph, morphine_RHP_graph, align = "h")
ggsave(filename = "Morphine_DR_plots.tiff", plot =Morphine_DR_plots, width = 10, height = 5)

## G*Power

### Alpha = 0.05, 0.01 SHP
ggplot(G_power, aes(x = Sample_size, y = logSHP, group = Groups)) +
  geom_point(aes(pch = factor(Groups), colour = factor(Groups)), size = 2, alpha = 0.75) +
  geom_line(aes(colour = factor(Groups))) +
  facet_grid(alpha ~ Treatment) +
  theme_bw() +
  scale_shape_manual(values = c(16, 15, 17, 18, 21), name = "Treatment\ngroups") +
  scale_x_continuous(limits = c(10, 80), breaks = seq(10, 80, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c("#1B9E77", "#404040", "#000000", 
                                "#E6AB02", "steelblue"), name = "Treatment\ngroups") +
  xlab("Total sample size") +
  ylab(expression("Power (1 - "*beta*" error probability)")) +
  labs(title = "G*Power analysis: Standard hot plate (log-transformed)",
       subtitle = expression("Rows determined by "*alpha*" value; Effect size: THC ("*omega^2*" = 0.67), morphine ("*omega^2*" = 0.71)")) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.x = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.y = element_text(size = 14, colour = "black", family = "Arial"),
        strip.background = element_rect(fill = alpha("seagreen", 0.75)),
        legend.text = element_text(size = 14, family = "Arial")) +
  geom_abline(slope = 0, intercept = 0.8, lty = 3, colour = "black", alpha = 0.75) +
  geom_abline(slope = 0, intercept = 0.9, lty = 1, colour = "black", alpha = 0.75)

### Combined G*Power plot

Gpower_plot_1 <- ggplot(Power_plot_2, aes(x = Sample_size, y = Value, group = Groups)) +
  geom_point(aes(pch = factor(Groups), colour = factor(Groups)), size = 2.5, alpha = 0.85) +
  geom_line(colour = "black", alpha = 0.5) +
  facet_grid(Protocol ~ Treatment) +
  theme_bw() +
  scale_shape_manual(values = c(16, 15, 17, 18, 21), name = "Treatment\ngroups") +
  scale_colour_grey(name = "Treatment\ngroups") +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ylab(expression("Power (1 - "*beta*" error probability)")) +
  xlab("Total sample size") +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.x = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0,
                                    family = "Arial"),
        strip.background.x = element_rect(fill = "transparent"),
        legend.text = element_text(size = 14, family = "Arial")) +
  geom_abline(slope = 0, intercept = 0.8, lty = 3, colour = "black", size = 0.5, alpha = 0.75) +
  geom_vline(data = mean.data, mapping = aes(xintercept = c(28, 25, 32, 23)),
             alpha = 0.75, size = .5) 

####
g1 <- ggplot_gtable(ggplot_build(Gpower_plot_1))
stripr1 <- which(grepl('strip-r', g1$layout$name))
fills1 <- alpha(c("#000000", "#404040"), 0.5)
k1 <- 1
for (i in stripr1) {
  j1 <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills1[k1]
  k1 <- k1+1
}
GPower_combined_plot1 <- plot_grid(g1)
ggsave(filename = "GPower_combined_plot2.tiff", plot = GPower_combined_plot1, width = 10, height = 7)

### Combined G*Power plot group size 5

Gpower_plot_2 <- ggplot(Power_plot_2, aes(x = Sample_size, y = Value, group = Groups)) +
  geom_point(aes(pch = factor(Groups), colour = factor(Groups)), size = 2.5, alpha = 0.85) +
  geom_point(data = subset(Power_plot_2, Groups == 5), aes(x = Sample_size, y = Value), 
             colour = "red", pch = 18, alpha = 0.5, size = 2.5) +
  geom_line(colour = "black", alpha = 0.5) +
  geom_line(data = subset(Power_plot_2, Groups == 5), aes(x = Sample_size, y = Value), 
            colour = "red", alpha = 0.5) +
  facet_grid(Protocol ~ Treatment) +
  theme_bw() +
  scale_shape_manual(values = c(16, 15, 17, 18, 21), name = "Treatment\ngroups") +
  scale_x_continuous(limits = c(10, 80), breaks = seq(10, 80, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_grey(name = "Treatment\ngroups") +
  xlab("Total sample size") +
  ylab(expression("Power (1 - "*beta*" error probability)")) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text = element_text(size = 14, colour = "black", family = "Arial"),
        strip.text.x = element_text(size = 14, colour = "white", family = "Arial"),
        strip.text.y = element_text(size = 14, colour = "black", angle = 0,
                                    family = "Arial"),
        strip.background = element_rect(fill = alpha("grey10", 0.75)),
        legend.text = element_text(size = 14, family = "Arial")) +
  geom_vline(data = mean.data, mapping = aes(xintercept = c(28, 25, 32, 23)),
             alpha = 0.75, size = .5) +
  geom_abline(slope = 0, intercept = 0.8, lty = 3, colour = "black", size = 0.5, alpha = 0.75)

#### Isolating 5 Group size 

g2 <- ggplot_gtable(ggplot_build(Gpower_plot_2))
stripr2 <- which(grepl('strip-r', g2$layout$name))
fills2 <- alpha(c("#000000", "#404040"), 0.5)
k2 <- 1
for (i in stripr2) {
  j2 <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j2]]$gp$fill <- fills2[k2]
  k2 <- k2+1
}

GPower_combined_plot2 <- plot_grid(g2)
ggsave(filename = "GPower_combined_plot2.tiff", plot = GPower_combined_plot2, width = 10, height = 7)


