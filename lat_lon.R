library(tidyverse)
library(bold)
library(getopt)

spec <- matrix(c(
  'taxon',    't', 2, 'character', 'BOLD search term',
  'output',   'o', 1, 'character', 'Output metadata csv'
), byrow = T, ncol = 5)

opt <- getopt(spec)

#getwd()
#setwd('~/OneDrive/Diversity/')

out <- bold_seqspec(taxon=opt$taxon)

print(paste(nrow(out), 'records found for', opt$taxon))

# Filter data into new dataframe with specified columns only
meta <- out %>% select(processid, bin_uri, markercode, family_taxID, family_name, 
                       subfamily_taxID, subfamily_name, genus_taxID, genus_name, 
                       species_taxID, species_name, lat, lon)
str(out)

# See how many NA in each column
sapply(meta, function(x) sum(is.na(x)))

# Split dataframe into list of dataframes, one for each 10 degrees of latitude
lats <- split(meta, cut(meta$lat, seq(-90, 90, 10), include.lowest=TRUE))

length(lats)
nrow(lats[[5]])

l = seq(-90, 90, 10)
l[1]

# New empty dataframe
df <- data.frame(Lat = character(), Rows = numeric(), BINs = numeric(), Species = numeric())

# Loop through list of dataframes, count entries, unique species and unique BINs
x = 1
for (lat in lats) {
  Lat = paste(l[x], 'to', l[x] + 10)
  Rows = nrow(lat)
  BINs = if (is.null(length(unique(lat$bin_uri)))) 0 else length(unique(lat$bin_uri))
  Species = length(unique(lat$species_taxID))
  df <- rbind(df, data.frame(Lat = Lat, Rows = Rows, BINs = BINs, Species = Species))
  x <- x + 1
  }

# Add new column for barcode:BIN ratio
df$Barcode_BIN_Ratio <- df$Rows/df$BINs

# Barcode: Species ratio
df$Barcode_Species_Ratio <- df$Rows/df$Species

write.csv(df, opt$output, row.names = FALSE)

