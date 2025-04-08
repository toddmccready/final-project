library(tidyverse)
library(gganimate)
library(magick)
library(latex2exp)


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

# Calculate fano factors
fano = data %>%
    inner_join(pars, by = join_by(sim_id)) %>%
    group_by(sim_id) %>%
    summarize(theory = (mu[1] + mu[1]^2 / iodisp[1]) / mu[1],
              empirical = var(count) / mean(count)) %>%
    mutate(diff = empirical/theory)

my_plot = ggplot(fano, aes(empirical)) +
    geom_density(fill = 'lightblue') +
    geom_rug() +
    xlim(0, 1.5) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    labs(title = TeX("Comparing the Theoretical v.s. Empirical Fano Factor"),
         y = "Density",
         x = TeX("$\\ \\frac{\\hat{F}_i}{F_i}$")) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12) 
    )
ggsave("figures/fano.png", my_plot,
       width = 7, height = 4)