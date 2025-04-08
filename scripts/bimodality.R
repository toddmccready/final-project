library(tidyverse)
library(mousetrap)
library(ggbeeswarm)


# Read in data
data = read.csv("data.csv")

# Remove beginning zeros to avoid non-stationary samples
data = data %>%
    group_by(sim_id) %>%
    mutate(row_number = row_number(),
           first_nonzero = min(row_number[count != 0], default = Inf)) %>%
    filter(row_number >= first_nonzero) %>%
    select(-row_number, -first_nonzero) %>%
    ungroup()

# Compute bimodality coefficient per-simulation
data = data %>%
    group_by(sim_id) %>%
    summarize(bc = bimodality_coefficient(count))

# Create distribution
my_plot = ggplot(data, aes(x = bc)) +
    geom_density(fill = 'lightblue', alpha = 0.5) +
    geom_rug() +
    xlim(0,1) +
    geom_vline(xintercept = 5/9, linetype = 'dashed') +
    theme_classic() +
    labs(x = "Bimodality Coefficient",
         y = "Density",
         title = "Distribution of Bimodality Coefficients Across Simulations") +
    theme(
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12) 
    )
ggsave("figures/bimodality.png", my_plot,
       width = 7, height = 1*3)
