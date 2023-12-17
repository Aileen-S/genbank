library(tidyverse)
library(bold)
library(getopt)

spec <- matrix(c(
  'taxon',  't', 2, 'character', 'BOLD search term',
  'input',  'i', 2, 'character', 'BOLD csv file (instead of --taxon)',
  'output', 'o', 1, 'character', 'Output metadata csv',
  'raw',    'r', 2, 'logical', 'Also save raw metadata before processing'
), byrow = T, ncol = 5)

opt <- getopt(spec)

#getwd()
#setwd('~/OneDrive/Diversity/')
out <- bold_seqspec(taxon='Eretes')


if ( !is.null(opt$taxon) ) {
  out <- bold_seqspec(taxon=opt$taxon)
} else {
  out <- read.csv(opt$input)
} 

if ( !is.null(opt$raw) ) {
  write.csv(out, paste('raw', opt$output, sep = '_'), row.names = FALSE)
}


print(paste(nrow(out), 'records found for', opt$taxon))

# Filter data into new dataframe with specified columns only
meta <- out %>% select(processid, bin_uri, markercode, family_taxID, family_name, 
                       subfamily_taxID, subfamily_name, genus_taxID, genus_name, 
                       species_taxID, species_name, lat, lon)
#str(out)

# See how many NA in each column
#sapply(meta, function(x) sum(is.na(x)))

# Split dataframe into list of dataframes, one for each 10 degrees of latitude
lats <- split(meta, cut(meta$lat, seq(-90, 90, 10), include.lowest=TRUE))

l = seq(-90, 90, 10)

# New empty dataframe
df <- data.frame(Lat = character(), IDs = numeric(), BINs = numeric(), Species = numeric())

# Loop through list of dataframes, count entries, unique species and unique BINs
x = 1
for (lat in lats) {
  Lat = paste(l[x], 'to', l[x] + 10)
  IDs = if (is.null(length(unique(lat$processid)))) 0 else length(unique(lat$processid))
  BINs = if (is.null(length(unique(lat$bin_uri)))) 0 else length(unique(lat$bin_uri))
  Species = length(unique(lat$species_taxID))
  df <- rbind(data.frame(Lat = Lat, IDs = IDs, BINs = BINs, Species = Species), df)
  x <- x + 1
  }

# Add new column for barcode:BIN ratio
df$Barcode_BIN_Ratio <- df$IDs/df$BINs

# Barcode: Species ratio
df$Barcode_Species_Ratio <- df$IDs/df$Species

df

write.csv(df, opt$output, row.names = FALSE)

