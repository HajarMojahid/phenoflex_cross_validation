
# load libraries

library(ggplot2)
library(tidyverse)
library(cowplot)


#### loading parameters #####

### Spain ####
#load param global opt  
param_old_es <- read.csv('old_fitter/08_summary_pars_ES.csv')
#organize 
param_old_es$fitter <- c('old')
#homogenize
param_old_es$Cultivar <- sub("^data_", "", param_old_es$Cultivar)


#load param enhanced global opt 
param_new_es <- read.csv ('new_fitter/ES_fit/param_sum_es.csv')
#organize 
param_new_es$fitter <- c('new')

### Germany ####

#load param global opt  
param_old_de <- read.csv('old_fitter/09_summary_pars_DE.csv')
#organize 
param_old_de$fitter <- c('old')
#homogenize
param_old_de$Cultivar <- sub("^data_", "", param_old_de$Cultivar)
param_old_de <- param_old_de %>%
  mutate(Cultivar = recode(Cultivar, 
                           RoterBoskoop = 'Roter Boskoop',
                           JamesGrieve = 'James Grieve',
                           GoldenDelicious = 'Golden Delicious',
                           CoxOrange = 'Cox Orange'))

#load param enhanced global opt 
param_new_de <- read.csv ('new_fitter/de_fit/param_sum_DE.csv')
#organize 
param_new_de$fitter <- c('new')


### Apple ####
#load param global opt 
param_old_apple <- read.csv('old_fitter/07_summary_pars_both.csv')
#organize 
param_old_apple$fitter <- c('old')

#homogenize
param_old_apple$Cultivar <- sub("^data_", "", param_old_apple$Cultivar)

param_old_apple <- param_old_apple %>%
  mutate(Cultivar = recode(Cultivar, 
                           RoterBoskoop = 'Roter Boskoop',
                           JamesGrieve = 'James Grieve',
                           GoldenDelicious = 'Golden Delicious',
                           CoxOrange = 'Cox Orange'))

#load param enhanced global opt 
param_new_apple <- read.csv ('new_fitter/de_es_fit/param_sum_de_es.csv')
#organize 
param_new_apple$fitter <- c('new')


##### common parameters ######

## Spain
common_param_old_es <- param_old_es %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

common_param_new_es <- param_new_es %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

##### Germany
common_param_old_de <- param_old_de %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

common_param_new_de <- param_new_de %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

##### Apple
common_param_old_apple <- param_old_apple %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

common_param_new_apple <- param_new_apple %>%
  filter(!(Par %in% c('yc', 'zc', 's1')))  %>%
  group_by(Par,Fit_res,fitter) %>% 
  summarise(Value = mean(Value)) # getting rid of repeated values

#put all parameters together 
#first add Fit  

common_param_es <- rbind(common_param_old_es,common_param_new_es)
common_param_es <- common_param_es %>%
  mutate( Fit ="Location specific, Spain")

common_param_de <- rbind(common_param_old_de, common_param_new_de)
common_param_de <- common_param_de %>%
  mutate(Fit ="Location specific, Germany") 

common_param_apple <- rbind(common_param_new_apple,common_param_old_apple)
common_param_apple <- common_param_apple %>%
  mutate(Fit = "Apple fit")

#bind
common_param_all <- rbind(common_param_es,
                          common_param_de,
                          common_param_apple)

#save
write.csv(common_param_all, "results/parameters/parameters_common_all_fits_runs_cultivars_fitters.csv",
          row.names = FALSE)

#read parameters dataframe
#common_param_all <- read.csv("results/parameters/parameters_common_all_fits_runs_cultivars_fitters.csv")

#prepare to plot 
# rename s to slope

# use quartiles 

common_param_summary_fitter <- common_param_all %>%
  mutate(Par = recode(Par, "s" = "slope")) %>% 
  group_by(Par, Fit, fitter) %>%
  summarise(median_value = median(Value),
            q25 = quantile(Value, 0.25),
            q75 = quantile(Value, 0.75),
            min_value = q25,
            max_value = q75)


# Define the order of parameters
param_order <- c("A0", "A1", "E0", "E1", "slope", "Tb", "Tc", "Tf", "Tu")

  
# Define the labeller function
param_labeller <- function(variable, value) {
  common_labels <- list(
    'A0' = expression(bold("A0")),
    'A1' = expression(bold("A1")),
    'E0' = expression(bold("E0")),
    'E1' = expression(bold("E1")),
    'slope' = expression(bold("slope")),
    'Tb' = expression(bold(T)[bold(b)]),
    'Tc' = expression(bold(T)[bold(c)]),
    'Tf' = expression(bold(T)[bold(f)]),
    'Tu' = expression(bold(T)[bold(u)])
  )
  return(common_labels[value])
}

# Create a factor with levels in the desired order
common_param_summary_fitter$Par <- factor(common_param_summary_fitter$Par,
                                          levels = param_order)

##### plot common parameters #####

p_common_param_summary_fitter <- ggplot(common_param_summary_fitter, 
                                        aes(x = median_value, y = Fit, color = Fit, shape = fitter)) +
  
  geom_point(size = 3, position = position_dodge(width = 0.8) ) +
  
  geom_errorbar(aes(xmin = min_value, xmax = max_value), position = position_dodge(width = 0.8),
                width = 0.2) +
  
  theme_bw(base_size = 14) +
  
  theme(
    
    # Facet formatting
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold.italic"),
    plot.margin = margin (1.5,1.5,4,1.5),
    panel.spacing.y = unit(0.2, "lines"), 
    panel.spacing.x = unit(0.2,"lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "gray20"),
    
    # Axis formatting
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    
    #legend formatting
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.01, "pt"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(t = -5, unit = "pt"),
    legend.box = "vertical",
    
  ) +
  
  facet_wrap(~ Par, scales = "free_x", ncol = 3, labeller = param_labeller) +
  
  # Colors and labels
  
  scale_color_manual(name = "Calibration approaches:",
                     values = c("#E64B35",  "#4DBBD5", "#BCD979" ),
                     labels = c("Species-specific","German location-specific","Spanish location-specific")) +
  
  scale_y_discrete(limits = rev) +  
  
  scale_shape_manual(values = c( "new" = 17, "old" = 16),
                     labels = c("E.G. optimization","G. optimization"),
                     name = "Fitting procedures:") +
  labs(x = "Parameter Value")

# Save the plot
ggsave("plots/parameters/common_parameters_fitter_quartiles.png",
       p_common_param_summary_fitter, 
       height = 20, width = 19, dpi = 500,
       units = 'cm', device = 'jpeg')



##### specific parameters #####

#load parameters 

### Spain 
specific_param_old_es <- param_old_es%>%
  filter (Par %in% c('yc', 'zc', 's1'))
specific_param_new_es <- param_new_es%>%
  filter (Par %in% c('yc', 'zc', 's1'))

### Germany
specific_param_old_de <- param_old_de%>%
  filter (Par %in% c('yc', 'zc', 's1'))
specific_param_new_de <- param_new_de%>%
  filter (Par %in% c('yc', 'zc', 's1'))

### Apple 
specific_param_old_apple <- param_old_apple%>%
  filter (Par %in% c('yc', 'zc', 's1'))
specific_param_new_apple <- param_new_apple%>%
  filter (Par %in% c('yc', 'zc', 's1'))

#put all parameters together 
#first add location and fit 

specific_param_es <- rbind(specific_param_old_es,specific_param_new_es)
specific_param_es <- specific_param_es %>%
  mutate(Location ="Spain", Fit ="Location specific")

specific_param_de <- rbind(specific_param_old_de, specific_param_new_de)
specific_param_de <- specific_param_de %>%
  mutate(Location ="Germany", Fit ="Location specific")

specific_param_apple <- rbind(specific_param_new_apple,specific_param_old_apple)
specific_param_apple <- specific_param_apple %>%
  mutate(Location = case_when(Cultivar %in% specific_param_es$Cultivar ~ "Spain",
                              TRUE ~ "Germany"), Fit = "Apple fit")
#bind
specific_param_all <- rbind(specific_param_es,
                            specific_param_de,
                            specific_param_apple) %>%
  mutate(Cultivar = recode(Cultivar, # some order 
                           RoterBoskoop = 'Roter Boskoop',
                           JamesGrieve = 'James Grieve',
                           GoldenDelicious = 'Golden Delicious',
                           CoxOrange = 'Cox Orange',
                           DelaRiega = 'De la Riega'))


#save
write.csv(specific_param_all, "results/parameters/parameters_specific_all_fits_runs_cultivars_fitters.csv",row.names = FALSE)

#read parameters
#specific_param_all <- read.csv("results/parameters/parameters_specific_all_fits_runs_cultivars_fitters.csv")


##prepare to plot 

# use quartiles 
specific_param_summary_fitter <- specific_param_all %>%
  group_by(Location, Cultivar, Par, Fit, fitter) %>%
  summarise(
    median_value = median(Value),
    q25 = quantile(Value, 0.25),
    q75 = quantile(Value, 0.75),
    min_value = q25,
    max_value = q75
  ) %>%
  mutate(
    Par = factor(Par, levels = c("yc", "zc", "s1"))
  )



###### plot specific par by facets German cultivars #####
#to solve labeling problem
## have common scale for both locations but different fo each facet 
# load library for setting 


library(ggh4x)

parameter_scales <- list(
  yc = scale_x_continuous(name = "Parameter value", limits = c(min(specific_param_summary_fitter$min_value[specific_param_summary_fitter$Par == "yc"]) - 0.5,
                                                               max(specific_param_summary_fitter$max_value[specific_param_summary_fitter$Par == "yc"]) + 0.5)),
  zc = scale_x_continuous(name = "Parameter value", limits = c(min(specific_param_summary_fitter$min_value[specific_param_summary_fitter$Par == "zc"]) - 0.5,
                                                               max(specific_param_summary_fitter$max_value[specific_param_summary_fitter$Par == "zc"]) + 0.5)),
  s1 = scale_x_continuous(name = "Parameter value", limits = c(min(specific_param_summary_fitter$min_value[specific_param_summary_fitter$Par == "s1"]) - 0.5,
                                                               max(specific_param_summary_fitter$max_value[specific_param_summary_fitter$Par == "s1"]) + 0.5))
)


  
# Define the labeller function
  param_labeller_sp <- function(variable, value) {
    if (variable == "Par") {
     sp_labels <- list(
      'yc' = expression(bold(y)[bold(c)]),
      'zc' = expression(bold(z)[bold(c)]),
      's1' = expression(bold(s)[bold(1)])
    )
    return(sp_labels[value])
  } else {
    # Default label for fitter
    return(value)
  }
}
  

german_parameter <- ggplot(filter(specific_param_summary_fitter, Location == "Germany"),
                           aes(x = median_value, y = Cultivar, color = Fit)) +
  # base layers
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(xmin = min_value, xmax = max_value),
                position = position_dodge(width = 1), 
                width = 0.2) +
  #facet by fitter 
  facet_grid(factor(fitter, levels = c("old","new"),
                    labels = c("G. optimization","E.G. optimization")) ~ 
             Par, scales = "free", labeller = param_labeller_sp) +
  
  facetted_pos_scales(x = parameter_scales) +
  
  theme_bw() +
  
  theme(legend.position = "none",
        axis.title.y =  element_text( size = 11, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold")
        ) +
  
  #set scale 
  scale_x_continuous(limits = parameter_scales) +
  
  # Colors and labels
  scale_color_manual(name = "Calibration approaches:",
                     values = c("Location specific" = "#4DBBD5", 
                                "Apple fit" = "#E64B35"),
                     labels = c( "Species-specific","Location-specific")) +
  
  labs(y = "Germany")

###### plot specific par  by facets Spanish cultivars #####

spanish_parameter <- ggplot(filter(specific_param_summary_fitter, Location == "Spain"),
                           aes(x = median_value, y = Cultivar, color = Fit)) +
  # base layers
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(xmin = min_value, xmax = max_value),
                position = position_dodge(width = 1), 
                width = 0.2) +
  #facet by fitter 
  facet_grid(factor(fitter, levels = c("old","new"), labels = c("G. optimization","E.G. optimization")) ~ Par,
             scales = "free") +
  #set scales
  facetted_pos_scales(x = parameter_scales) +
  
  theme_bw() +
  
  theme(legend.position = "bottom",
        legend.title = element_text(),
        legend.margin = margin(t = -10, unit = "pt"),
        axis.title.y =  element_text( size = 11, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text( size = 11, face = "bold")
        
  ) +
  
  
  # Colors and labels
  scale_color_manual(name = "Calibration approaches:",
                     values = c("Location specific" = "#4DBBD5", 
                                "Apple fit" = "#E64B35"),
                     labels = c("Species-specific","Location-specific"))+
  
  labs(y = "Spain")  

##### use cowplot to put the plots together #####


combined_plot <- plot_grid(german_parameter, spanish_parameter,
          ncol = 1, align = "v",rel_heights = c (1,2))

# add cultivar on y-axis 

final_plot <- ggdraw(combined_plot) +
  draw_label("Cultivar", x = 0.08, y = 0.97, hjust = 0, vjust = 0, fontface = "bold", size = 10) 


## save 
ggsave("plots/parameters/specific_parameters_dpi.png", 
       height = 24, width = 19, dpi = 500,
       units = 'cm', device = 'jpeg')
