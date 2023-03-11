library(faux)
library(tidyverse)


cmat <- c(1, .3,
          .3, 1)

df <- rnorm_multi(n = 200, vars = 2, mu = 0, sd = c(1,2), cmat, 
                   varnames = c("g", "s")) %>%
  mutate("group" = rep(1:4, each = 50) ** 2,
         "site" = rep(1:4, 50), 
         "x" = group + rnorm(200,0,0.5),
         "y" = rnorm(200, x/2, sqrt(abs(x))) + g + s,
         group = factor(group),
         site = factor(site))

# Plot phenotypes

df %>%
  ggplot(aes(x = x, y = y, col = site)) +
  geom_point()


df %>%
  mutate(res = resid(lm(y~ x, data = .))) %>%
  ggplot(aes(x = x, y = res, col = site)) +
  geom_point()
