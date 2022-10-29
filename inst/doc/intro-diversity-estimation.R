## -----------------------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("adw96/breakaway")
library(breakaway)

## -----------------------------------------------------------------------------
# install.packages("magrittr")
library(magrittr)

## -----------------------------------------------------------------------------
# install.packages("openxlsx")
# pasolli_et_al <- openxlsx::read.xlsx("https://ars.els-cdn.com/content/image/1-s2.0-S0092867419300017-mmc4.xlsx", sheet = 2)
data("pasolli_et_al")

## -----------------------------------------------------------------------------
ft <- pasolli_et_al$`#.Samples` %>% make_frequency_count_table
ft %>% head(10)

## -----------------------------------------------------------------------------
ft %>% tail(10)

## -----------------------------------------------------------------------------
ft %>% sample_richness

## -----------------------------------------------------------------------------
estimated_richness <- breakaway(ft)
estimated_richness

## -----------------------------------------------------------------------------
estimated_richness$model

## -----------------------------------------------------------------------------
plot(estimated_richness)

## -----------------------------------------------------------------------------
estimated_richness_10 <- breakaway(ft, cutoff = 10)
estimated_richness_10

## -----------------------------------------------------------------------------
estimated_richness_10 %$% model

## -----------------------------------------------------------------------------
estimated_richness_kemp_10 <- kemp(ft, cutoff = 10)
estimated_richness_kemp_10

## -----------------------------------------------------------------------------
estimated_richness_kemp_10 %$% plot

