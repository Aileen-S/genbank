library('getopt')
library(phytools)

spec <- matrix(c(
  'input',  'i', 2, 'character', 'Input tree',
  'output', 'o', 2, 'character', 'Output binary tree'
), byrow = T, ncol = 5)
opt <- getopt(spec)

tree <- read.tree(opt$input)

# Check if the tree is binary
if (is.binary(tree)) {
  cat(opt$input, 'is binary\n')
} else {
  tree <- multi2di(tree)
  write.tree(tree, file=opt$output)
  cat(opt$input, 'is not binary\n')
  cat('Binary tree written to', opt$output, '\n')
}
