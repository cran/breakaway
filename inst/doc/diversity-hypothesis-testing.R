## ---- include=FALSE-----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE) 

## ---- include=FALSE-----------------------------------------------------------
library(tidyverse)
library(breakaway)
library(magrittr)
library(phyloseq)

## -----------------------------------------------------------------------------
library(corncob)
data("soil_phylo")
soil_phylo %>% sample_data %>% head

## -----------------------------------------------------------------------------
subset_soil <- soil_phylo %>%
  subset_samples(Amdmt == 1) %>% # only biochar
  subset_samples(Day %in% c(0, 2))  # only Days 0 and 82

## -----------------------------------------------------------------------------
richness_soil <- subset_soil %>% breakaway
plot(richness_soil, physeq=subset_soil, color="Day", shape = "ID")

## -----------------------------------------------------------------------------
summary(richness_soil) %>% as_tibble

## -----------------------------------------------------------------------------
meta <- subset_soil %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = subset_soil %>% sample_names )

## -----------------------------------------------------------------------------
combined_richness <- meta %>%
  left_join(summary(richness_soil),
            by = "sample_names")
# Old way (still works)
bt_day_fixed <- betta(chats = combined_richness$estimate,
                      ses = combined_richness$error,
                      X = model.matrix(~Day, data = combined_richness))
# Fancy new way -- thanks to Sarah Teichman for implementing!
bt_day_fixed <- betta(formula = estimate ~ Day, 
                      ses = error, data = combined_richness)
bt_day_fixed$table

## -----------------------------------------------------------------------------
# Old way (still works)
bt_day_fixed_id_random <- betta_random(chats = combined_richness$estimate,
                                       ses = combined_richness$error,
                                       X = model.matrix(~Day, data = combined_richness),
                                       groups=combined_richness$ID)
# Fancy new way 
bt_day_fixed_id_random <-
  betta_random(formula = estimate ~ Day | ID, 
               ses = error,  data = combined_richness)
bt_day_fixed_id_random$table

## ---- include = FALSE---------------------------------------------------------
# remotes::install_github("adw96/DivNet")
# library(DivNet)

## -----------------------------------------------------------------------------
soil_phylum <- subset_soil %>%
  tax_glom(taxrank="Phylum")

## ---- include=TRUE, cache=TRUE, message=FALSE, results="hide"-----------------
# the following line contains the DivNet command to fit the model we will use
# dv <- DivNet::divnet(soil_phylum, X = NULL)

# to allow readers to build this vignette without installing DivNet
# we've saved the fitted model, which we'll load directly:

data("dv")

## -----------------------------------------------------------------------------
dv

## -----------------------------------------------------------------------------
combined_shannon <- meta %>%
  dplyr::left_join(dv$shannon %>% summary,
            by = "sample_names")
combined_shannon

## -----------------------------------------------------------------------------
bt_day_fixed_id_random <- betta_random(formula = estimate ~ Day | ID, 
               ses = error,  data = combined_shannon)
bt_day_fixed_id_random$table

## -----------------------------------------------------------------------------

betta_lincom(fitted_betta = bt_day_fixed_id_random,
             linear_com = c(1,1),
             signif_cutoff = 0.05)



## -----------------------------------------------------------------------------

subset_soil_days_1_2 <- soil_phylo %>%
  subset_samples(Amdmt == 1) %>% # only biochar
  subset_samples(Day %in% c(0, 1, 2))  # only Days 0 and 82

## -----------------------------------------------------------------------------
meta_days_1_2 <- subset_soil_days_1_2 %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = subset_soil_days_1_2 %>% sample_names )

soil_phylum_days_1_2 <- subset_soil_days_1_2 %>%
  tax_glom(taxrank="Phylum")

## -----------------------------------------------------------------------------
# we again have stored the output of the following DivNet call
# dv_days_1_2 <- DivNet::divnet(soil_phylum_days_1_2, X = NULL)

# and we load the fitted model directly
data("dv_days_1_2")

combined_shannon_days_1_2 <- meta_days_1_2 %>%
  dplyr::left_join(dv_days_1_2$shannon %>% summary,
            by = "sample_names")
combined_shannon_days_1_2


## -----------------------------------------------------------------------------
bt_day_1_2_fixed_id_random <- betta_random(formula = estimate ~ Day | ID, 
               ses = error,  data = combined_shannon_days_1_2)
bt_day_1_2_fixed_id_random$table

## -----------------------------------------------------------------------------
set.seed(345)
submodel_test <- test_submodel(bt_day_1_2_fixed_id_random,
                          submodel_formula = estimate~1,
                          method = "bootstrap",
                          nboot = 100)

submodel_test$pval


