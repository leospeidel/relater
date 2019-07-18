# relater

<!-- badges: start -->
<!-- badges: end -->

R package for parsing output files generated by Relate, as well as simple manipulations of them, including an implementation of the polygenic selection test.

Please note that this package is still under development, although its main functionalities should be working.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("leospeidel/relater")
```

## Example

``` r
library(relater)

################################################
# Example analysis for polygenic selection test:

# read files
mut         <- read.mut(system.file("extdata/example.mut.gz", package = "relater"), CHR = 1)
sele        <- read.sele(system.file("extdata/example.sele.gz", package = "relater"), CHR = 1)
freq        <- read.freq(system.file("extdata/example.freq.gz", package = "relater"), CHR = 1)
qual        <- read.qual(system.file("extdata/example.qual.gz", package = "relater"), CHR = 1)


# Obtain allele_ages data table
allele_ages <- get.allele_ages(mut, freq, sele)
allele_ages <- filter.allele_ages(allele_ages, qual)

######### Polygenic selection test #########

# Initialise
quant       <- PolyTest_Init(allele_ages)

# Make a fake polygenic trait
df          <- allele_ages[!is.na(allele_ages$pvalue),]
df          <- df[df$pvalue < 0,]
df          <- df[sample(1:nrow(df), 50, replace = FALSE),]

# Run polygenic test
PolyTest(df$DAF, df$pvalue, quant$quantiles, quant$allele_ages_quantiles)

################################################
# Parsing .coal files

coal <- read.coal(system.file("extdata/example.coal.gz", package = "relater"))
head(coal)

################################################
# Parsing *_avg.rate files

avg_mutrate <- read.avg_rate(system.file("extdata/example_avg.rate.gz", package = "relater"))
head(avg_mutrate)

################################################
# Parsing haps/sample files

haps   <- read.haps(system.file("extdata/example.haps.gz", package = "relater"))
head(haps)

sample <- read.sample(system.file("extdata/example.sample.gz", package = "relater"))
head(sample)

```
