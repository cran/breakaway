## -----------------------------------------------------------------------------
library(breakaway)
data(toy_otu_table)
## For historical reasons we're going to rename it:
otu_data <- toy_otu_table

## -----------------------------------------------------------------------------
frequencytablelist <- build_frequency_count_tables(otu_data)
head(frequencytablelist[[63]])

## -----------------------------------------------------------------------------
breakaway(frequencytablelist[[1]])

## -----------------------------------------------------------------------------
plot(breakaway(frequencytablelist[[1]]))

## -----------------------------------------------------------------------------
breakaway(frequencytablelist[[2]])

## -----------------------------------------------------------------------------
breakaway_nof1(frequencytablelist[[2]][-1,])

