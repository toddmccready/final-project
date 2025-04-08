library(tidyverse)
library(gganimate)
library(magick)


# Functions ----
ll <- function(params, count) {
  log_mu <- params[1]
  log_iodisp <- params[2]

  mu <- exp(log_mu)
  iodisp <- exp(log_iodisp)

  log_probs <- dnbinom(count, mu = mu, size = iodisp, log = TRUE)
  return(-sum(log_probs))
}


# Main ----
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

# Compute best negative binomial parameters
pars = lapply(unique(data$sim_id), function(i) {
    cur_data = data %>% filter(sim_id == i)
    
    pars = optim(
        c(0,0),
        ll,
        count = cur_data$count
    )$par |> exp()
    
    tibble(
        sim_id = i,
        mu = pars[1],
        iodisp = pars[2]
    )
})
pars = do.call(rbind, pars)

# Include parameters
data = data %>%
    inner_join(pars, by = join_by(sim_id))

# Subset simulations for plotting
data = data %>%
    filter(sim_id < 30)

# Get the range of count values to plot the PMF over
count_range = min(data$count):max(data$count)

# Create a dataframe to store the PMF for each sim_id
pmf_data = data %>% 
    select(-count) %>% 
    unique() %>%
    group_by(sim_id) %>%
    crossing(tibble(
        count = count_range
    )) %>%
    mutate(prob = dnbinom(count, mu = mu, size = iodisp))

# Create animation
animated_plot = ggplot(data, aes(x = count, color = "True Distribution")) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.8) +
    geom_line(data = pmf_data, aes(x = count, y = prob, color = "Negative Binomial Approximation", group = sim_id)) +
    labs(title = "Negative Binomial Approximation to the True Distribution",
        x = "Count",
        y = "Probability Mass/Observed Frequency",
        color = "") +
    scale_color_manual(values = c("Negative Binomial Approximation" = "red")) +
    theme_classic() +
    transition_manual(sim_id) +
    theme(
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12)
    )


# Save animation
animation = animate(animated_plot, 
        nframes = length(unique(data$sim_id)), 
        fps = 2, 
        res = 300,
        height = 160*7.5, width = 160*16,
        renderer = magick_renderer())
anim_save("figures/approxs.gif", animation)
