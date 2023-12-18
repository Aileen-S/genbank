library(tidyverse)
library(bold)
library(getopt)

uk_total = data.frame(Family = c("Aderidae", "Alexiidae", "Anthicidae", "Anthribidae", "Apionidae", "Attelabidae", "Biphyllidae", "Bostrichidae", "Bothrideridae", "Buprestidae", "Byrrhidae", "Byturidae", "Cantharidae", "Carabidae", "Cerambycidae", "Cerylonidae", "Chrysomelidae", "Ciidae", "Clambidae", "Cleridae", "Coccinellidae", "Colydiidae", "Corylophidae", "Cryptophagidae", "Cucujidae", "Curculionidae", "Dascillidae", "Dasytidae", "Dermestidae", "Derodontidae", "Drilidae", "Dryophthoridae", "Dryopidae", "Dytiscidae", "Elateridae", "Elmidae", "Endomychidae", "Erirhinidae", "Erotylidae", "Eucinetidae", "Eucnemidae", "Georissidae", "Geotrupidae", "Gyrinidae", "Haliplidae", "Helophoridae", "Heteroceridae", "Histeridae", "Hydraenidae", "Hydrochidae", "Hydrophilidae", "Hygrobiidae", "Kateretidae", "Laemophloeidae", "Lampyridae", "Latridiidae", "Leiodidae", "Limnichidae", "Lucanidae", "Lycidae", "Lymexylidae", "Malachiidae", "Megalopodidae", "Melandryidae", "Meloidae", "Monotomidae", "Mordellidae", "Mycetophagidae", "Mycteridae", "Nanophyidae", "Nemonychidae", "Nitidulidae", "Noteridae", "Oedemeridae", "Orsodacnidae", "Phalacridae", "Phloiophilidae", "Platypodidae", "Psephenidae", "Ptiliidae", "Ptilodactylidae", "Ptinidae", "Pyrochroidae", "Pythidae", "Raymondionymidae", "Rhynchitidae", "Ripiphoridae", "Salpingidae", "Scarabaeidae", "Scirtidae", "Scraptiidae", "Silphidae", "Silvanidae", "Spercheidae", "Sphaeritidae", "Sphaeriusidae", "Sphindidae", "Staphylinidae", "Tenebrionidae", "Tetratomidae", "Throscidae", "Trogidae", "Trogossitidae"),
                      UK_Species = c(3, 1, 13, 9, 90, 2, 2, 10, 5, 18, 13, 2, 42, 374, 69, 5, 286, 22, 10, 15, 53, 12, 11, 105, 2, 493, 1, 9, 40, 1, 1, 4, 9, 120, 73, 12, 8, 13, 8, 1, 7, 1, 8, 12, 19, 20, 8, 52, 34, 7, 72, 1, 9, 11, 3, 56, 95, 1, 4, 3, 2, 17, 3, 17, 11, 24, 17, 15, 1, 2, 1, 91, 2, 10, 2, 15, 1, 2, 1, 75, 1, 57, 3, 1, 1, 19, 1, 11, 88, 20, 17, 21, 12, 1, 1, 1, 2, 1130, 47, 4, 5, 3, 5))

# Function to retrieve the number of species for a specific family
getSpeciesCount <- function(family) {
  subset_data <- subset(uk_total, Family == family)
  if (nrow(subset_data) > 0) {
    return(subset_data$UK_Species)
  } else {
    return("Family not found")
  }
}


spec <- matrix(c(
  'taxon',   't', 1, 'character', 'Coleoptera family',
  'input',   'i', 2, 'character', 'input BOLD csv file (optional, otherwise will search BOLD)',
  'output',  'o', 1, 'character', 'Output file',
  'raw',     'r', 2, 'logical',   'Also save raw metadata before processing',
  'country', 'c', 2, 'character', 'Limit output to specified country'
), byrow = T, ncol = 5)

opt <- getopt(spec)


#out <- bold_seqspec(taxon='Longitarsus')
#out <- read.csv('raw_Dytiscidae.csv')


if ( !is.null(opt$input) ) {
  out <- read.csv(opt$input)
} else {
  out <- bold_seqspec(taxon=opt$taxon)
}


if ( !is.null(opt$raw) ) {
  write.csv(out, paste('raw', opt$output, sep = '_'), row.names = FALSE)
}

print(paste(nrow(out), 'records found'))

# Remove records without BIN
meta <- out[!is.na(out$bin_uri),]
print(paste(nrow(meta), 'records with BIN'))

# Remove records without UK
meta <- meta %>% filter(country == opt$country)
print(paste(nrow(meta), 'records from UK'))

barcodes <- nrow(meta)
print(paste(barcodes, 'barcodes'))

BINs <- length(unique(meta$bin_uri))
print(paste(BINs, 'BINs'))

bar_bin <- barcodes/BINs
print(paste('Barcode:BIN ratio', bar_bin))

spec <- getSpeciesCount(opt$taxon)
coverage <- BINs/spec
print(paste(coverage, 'of UK species represented (approx)'))

sink(opt$output)
cat(paste(barcodes, BINs, bar_bin, spec, coverage, sep = ','))
sink()
print(paste('Stats written to', opt$output))

