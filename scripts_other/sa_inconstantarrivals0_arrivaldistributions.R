library(tidyverse)
library(gridExtra)

# plotting theme
my_theme <- function(){
  theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18, face = "bold"),
          plot.title=element_text(size=20, face = "bold", hjust = 0.5))
}

# Values from which to extract from the binomial distributions for B and C
vals <- c(.0625, .1875, .3125, .4375, .5625, .6875, .8125, .9375)

# Unimodal distribution 
dist_a <- c(rep(1/6, 4), rep(1/2,2), rep(1/6,2))
plot(vals*4, dist_a, ylim = c(0, 0.5))
a_df <- data.frame(x = vals*4, y = dist_a)
p1 <- ggplot(data = a_df, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  xlab("Hour") +
  ylab("Arrival density") +
  ggtitle("Density A") +
  my_theme()

# Unimodal distribution, skewed right
dist_b <- sapply(vals, function(x) dbeta(x, 2, 3)/4)
dist_b[3] <- dist_b[3] - sum(dist_b) + 2  # normalize values
plot(vals*4, dist_b, ylim = c(0, 0.5))
b_df <- data.frame(x = vals*4, y = dist_b)
p2 <- ggplot(data = b_df, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  xlab("Hour") +
  ylab("Arrival density") +
  ggtitle("Density B") +
  my_theme()

# Unimodal distribution, skewed left
dist_c <- sapply(vals, function(x) dbeta(x, 3, 2)/4)
dist_c[6] <- dist_c[6] - sum(dist_c) + 2  # normalize values
plot(vals*4, dist_c, ylim = c(0, 0.5))
c_df <- data.frame(x = vals*4, y = dist_c)
p3 <- ggplot(data = c_df, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  xlab("Hour") +
  ylab("Arrival density") +
  ggtitle("Density C") +
  my_theme()

# Bimodal distribution
dist_d <- c(0.15, 0.4, 0.25, 0.2, 0.2, 0.25, 0.4, 0.15)
plot(vals*4, dist_d, ylim = c(0, 0.5))
d_df <- data.frame(x = vals*4, y = dist_d)
p4 <- ggplot(data = d_df, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  xlab("Hour") +
  ylab("Arrival density") +
  ggtitle("Density D") +
  my_theme()

grid.arrange(p1, p2, p3, p4, ncol = 2)

# Convert values to /minute (instead of /hour)
a_min <- dist_a/60
b_min <- dist_b/60
c_min <- dist_c/60
d_min <- dist_d/60

# Save values 
dist_list <- list(dist_a = a_min, dist_b = b_min, dist_c = c_min, dist_d = d_min)
write_rds(dist_list, "data/arrivalratedistributions.rds")
