library(tidyverse)
library(lubridate)
library(gridExtra)

setwd("~/Optimization/")

dt <- read.csv("data/queue_data_2017.csv")

table(dt$LUGAR, dt$FECHA)

# Convert to time format
dt$time_arrival <- parse_date_time(dt$HORA_LLEGADA, "%H:%M")
dt$time_arrival[297:456] <- parse_date_time(dt$HORA_LLEGADA[297:456], "%H:%M:%S")
dt$time_start[297:456] <- parse_date_time(dt$INICIO_ATENCION[297:456], "%H:%M:%S")
dt$time_finish[297:456] <- parse_date_time(dt$FIN_ATENCION[297:456], "%H:%M:%S")
dt$time_depart[297:456] <- parse_date_time(dt$HORA_SALIDA[297:456], "%H:%M:%S")

# 1 failed to parse
which(is.na(dt$time_finish))
dt$FIN_ATENCION[405]

# Calculate service time and system time
dt$service_time <- dt$time_finish - dt$time_start

sum(is.na(dt$time))
which(is.na(dt$time))
dt$HORA_LLEGADA[is.na(dt$time)]

# Remove duplicates 
dt <- dt %>%
  select(- ID_numero) %>%
  distinct()

# Subset data into observation periods
dt1 <- filter(dt, LUGAR == "ESQUINA DE LA 25")
dt2 <- filter(dt, LUGAR == "ÓVALO HUNTER")
dt3 <- filter(dt, FECHA == "12/8/17")
dt4 <- filter(dt, FECHA == "13/08/17")

# plotting theme
my_theme <- function(){
  theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18, face = "bold"),
          plot.title=element_text(size=20, face = "bold", hjust = 0.5))
}

# Make arrival time histograms
# Round down to nearest 30 minute interval
dt1$time_floor <- floor_date(dt1$time_arrival, "30 mins")
p1 <- ggplot(data = dt1, aes(x = time_floor)) +
  geom_bar() +
  xlab("Time") +
  ylab("Number of arrivals") +
  ggtitle("Esquina de la 25, Day 1") +
  my_theme()

dt2$time_floor <- floor_date(dt2$time_arrival, "30 mins")
p2 <- ggplot(data = dt2, aes(x = time_floor)) +
  geom_bar() +
  xlab("Time") +
  ylab("Number of arrivals") +
  ggtitle("Óvalo Hunter, Day 2") +
  my_theme()

dt3$time_floor <- floor_date(dt3$time_arrival, "30 mins")
p3 <- ggplot(data = dt3, aes(x = time_floor)) +
  geom_bar() +
  xlab("Time") +
  ylab("Number of arrivals") +
  ggtitle("Microred San Martín, Day 2") +
  my_theme()

dt4$time_floor <- floor_date(dt4$time_arrival, "30 mins")
p4 <- ggplot(data = dt4, aes(x = time_floor)) +
  geom_bar() +
  xlab("Time") +
  ylab("Number of arrivals") +
  ggtitle("Microred San Martín, Day 3") +
  my_theme()

grid.arrange(p1, p2, p3, p4, ncol = 2)
