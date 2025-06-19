
library(dplyr)
library(ggplot2)
library(maps)

# locations
locations <- data.frame(
  name = c("CKA", "SERIDA"),
  lat = c(50.6254619, 43.4785112),
  lng = c(6.988372, -5.4389411)
)

# Base world map (Mediterranean/Europe focus)
world_BG <- borders("world", 
                    size = 0.2, 
                    fill = "lightgoldenrod2", 
                    colour = "black") 

# Create the plot
world <- ggplot() + 
  world_BG +
  geom_point(
    data = locations,
    aes(x = lng, y = lat),
    color = "black",    
    fill = "grey0",     
    size = 3,           
    shape = 21          
  ) +
  coord_sf(xlim = c(-10, 40), ylim = c(20, 70)) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.grid = element_blank())

# save 
ggsave(
  "plots/map.png", 
  plot = world,
  width = 190,      
  height = 200,     
  units = "mm",
  dpi = 500,
  device = "png"
)

# weather charts for each location 

# CKA 
clim_cka$Month <- factor(month.abb[clim_cka$Month], levels = month.abb)

p_cka <- ggplot(clim_cka, aes(x = Month)) +
  geom_bar(aes(y = Prec, fill = "Precipitation"), stat = "identity", alpha = 0.7) +
  geom_line(
    aes(y = Tmean * (max(Prec, na.rm = TRUE) / max(Tmean, na.rm = TRUE)), 
        color = "Temperature", group = 1) , linewidth = 0.5
  ) +
  geom_point(
    aes(y = Tmean * (max(Prec, na.rm = TRUE) / max(Tmean, na.rm = TRUE)), 
        color = "Temperature"), size = 0.25
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Precipitation" = "#4DBBD5"),
    labels = c("Precipitation (mm)")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "#E64B35"),
    labels = c("Temperature (°C)")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "#E64B35"),
    labels = c("Temperature (°C)")
  ) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    limits = c(0, 80),
    breaks = seq(0, 70, by = 10),  
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / (80 / 25),  
      name = "Temperature (°C)",
      breaks = seq(0, 20, by = 5)  
    )
  ) +
  labs(title = "CKA, Germany") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray20", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9, face = "bold"),
    plot.margin = margin(1.5, 1.5, 1.5, 2),
    plot.background = element_rect(colour = "white", fill = "white", size = 0.2),
    plot.title = element_text(size = 9, 
                              face = "bold",
                              hjust = 0.5,  
                              margin = margin(b = 5)),
    plot.title.position = "panel"
  ) +
  guides(color = "none", fill = "none")

# SERIDA 
clim_serida$Month <- factor(month.abb[clim_serida$Month], levels = month.abb)

p_serida <- ggplot(clim_serida, aes(x = Month)) +
  geom_bar(aes(y = Prec, fill = "Precipitation"), stat = "identity", alpha = 0.7) +
  geom_line(
    aes(y = Tmean * (max(Prec, na.rm = TRUE) / max(Tmean, na.rm = TRUE)), 
        color = "Temperature", group = 1) , linewidth = 0.5
  ) +
  geom_point(
    aes(y = Tmean * (max(Prec, na.rm = TRUE) / max(Tmean, na.rm = TRUE)), 
        color = "Temperature"), size = 0.25
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Precipitation" = "#4DBBD5"),
    labels = c("Precipitation (mm)")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "#E64B35"),
    labels = c("Temperature (°C)")
  ) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    limits = c(0, 80),
    breaks = seq(0, 70, by = 10),  
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / (80 / 25),  
      name = "Temperature (°C)",
      breaks = seq(0, 20, by = 5)  
    )
  ) +
  labs(title =  "SERIDA, Spain") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray20", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9, face = "bold"),
    plot.margin = margin(1.5, 1.5, 1.5, 2),
    plot.background = element_rect(colour = "white", fill = "white", size = 0.2),
    plot.title = element_text(size = 9, 
                              face = "bold",
                              hjust = 0.5,  
                              margin = margin(b = 5)),
    plot.title.position = "panel"
  ) +
  
  guides(color = "none", fill = "none")

# save 
ggsave(
  "plots/p_cka.png",
  plot = p_cka,
  width = 90, 
  height = 60, 
  units = "mm",
  dpi = 500,
  device = "png"
)
ggsave(
  "plots/p_serida.png",
  plot = p_serida,
  width = 90, 
  height = 60, 
  units = "mm",
  dpi = 500,
  device = "png"
)
