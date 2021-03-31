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

BPS_G_power_plots_copy <- read_excel("~/Desktop/Research/Cannevert /British Pharmacological Society/BPS abstract 2018/BPS 2018 presentation/Excel files/BPS G_power plots copy.xlsx")
GPower <- read_excel("~/Desktop/Research/Hot plate assay paper/Data analysis/GPower.xlsx")

## Processing files

### Control subpopulations

### SHP QQ-density plot

### RHP QQ-density plot

### Dose-response studies: Reference compounds

#### Standard hot plate

##### Pure THC

##### Morphine

#### Ramped hot plate

##### Pure THC

##### Morphine

### ANOVA

#### SHP

##### THC

##### Morphine

#### RHP

##### THC

##### Morphine

### G*Power

#### Alpha = 0.05, 0.01 SHP

GPower$treatment <- factor(GPower$treatment, ordered = T, levels = c("Morphine", "THC"))

#### Combined G*Power plot

x2_1 <- rep(1, 477)
x2_2 <- rep(2, 477)
x2_3 <- rep(3, 477)
x2_4 <- rep(4, 477)

binded_x2_ <- cbind(x2_1,
                   x2_2,
                   x2_3,
                   x2_4) 

binded2_df <- as.data.frame(binded_x2_)
melted2_df <- melt(binded2_df) %>% 
  dplyr::select(-variable)
melted2_df$value <- factor(melted2_df$value, levels = c(1, 2, 3, 4))

Combined_GPower <- GPower %>%
  gather(key = "Protocol", value = "Value", -sample_size, -groups, -treatment) %>%
  mutate(Protocol = gsub("SHP", "Standard\nhot plate", Protocol),
         Protocol = gsub("RHP", "Ramped\nhot plate", Protocol)) %>%
  dplyr::rename(Sample_size = sample_size,
                Treatment = treatment,
                Groups = groups)
  

Combined_GPower$Protocol <- factor(Combined_GPower$Protocol,
                                       levels = c("Standard\nhot plate", "Ramped\nhot plate"))

combined_melt2 <- cbind(Combined_GPower, melted2_df)
GPower_plot <- as.data.frame(combined_melt2)
GPower_plot_2 <- GPower_plot %>%
  mutate(line_SHP_THC = ifelse(value == 1, 32, NA),
         line_RHP_THC = ifelse(value == 2, 28, NA),
         line_SHP_morphine = ifelse(value == 3, 23, NA),
         line_RHP_morphine = ifelse(value == 4, 25, NA)) %>%
  gather(key = "factored_lines", value = "FLine",
         -Sample_size, -Groups, -Treatment, -Protocol, -Value, -value)

GPower_plot_2$FLine <- factor(GPower_plot_2$FLine, levels = c("32", "28", "23", "25"))

vline2_df <- data.frame(z = levels(GPower_plot_2$Treatment),
                       vl = c(32, 28, 23, 25))

##### Grouping data in order to allow for multiple lines specific to each panel

mean.data2 <- aggregate(x = GPower_plot_2$Value, # use the y values
                       by = GPower_plot_2[c("Protocol", "Treatment")], # group by Group and then by Subgroup
                       FUN = function(x) {
                         signif(mean(x), 4) # calculate the mean keep 4 significant numbers
                       }
)

colnames(mean.data2) <- c("Protocol", "Treatment", "Average")

# Figures

## Control subpopulations

### Image: Control subpopulation

## SHP QQ-density plot

### Untransformed

### Log-transformed

### Image: SHP QQ-density plot

## RHP QQ-density plot

### Untransformed 

### Log-transformed

### Image: RHP QQ-density plot

## Dose response curves: reference compounds

### Standard hot plate

#### Pure THC

#### Morphine

### Ramped hot plate

#### Pure THC

#### Moprhine

### Images

#### Pure THC

#### Morphine

## G*Power

### Alpha = 0.05, 0.01 SHP

### Combined G*Power plot

GPower_plot_4 <- ggplot(GPower_plot_2, aes(x = Sample_size, y = Value, group = Groups)) +
  geom_point(aes(pch = factor(Groups), colour = factor(Groups)), size = 2.5, alpha = 0.85) +
  geom_line(colour = "#202020", alpha = 0.5) +
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
        strip.text.x = element_text(size = 14, colour = "#404040", family = "Arial"),
        strip.text.y = element_text(size = 14, colour = "white", angle = 0,
                                    family = "Arial"),
        strip.background.x = element_rect(fill = "transparent"),
        legend.text = element_text(size = 14, family = "Arial")) +
  geom_abline(slope = 0, intercept = 0.8, lty = 3, colour = "black", size = 0.5, alpha = 0.75) 

####
gplot1 <- ggplot_gtable(ggplot_build(GPower_plot_4))
stripr1 <- which(grepl('strip-r', gplot1$layout$name))
fills1 <- alpha(c("#000000", "#404040"), 0.5)
k1 <- 1
for (i in stripr1) {
  j1 <- which(grepl('rect', gplot1$grobs[[i]]$grobs[[1]]$childrenOrder))
  gplot1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills1[k1]
  k1 <- k1+1
}
GPower_color_plot <- plot_grid(gplot1)
ggsave(filename = "GPower_combined_plot1.tiff", plot = GPower_color_plot, width = 10, height = 7)

### Combined G*Power plot group size 5


