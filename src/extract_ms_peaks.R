library(BiocParallel)
register(SerialParam())
library(msdata)
f <- c("~/data/20210701_10043/test/210705_AR_06_V_03_2_01.mzML.gz")

## Read the data as an MSnExp
msd <- readMSData(f, msLevel = 1)

## Extract chromatograms for a MS data slices defined by retention time
## and mz ranges.
rtr <- rbind(c(0, 1000), c(0, 1000))
mzr <- rbind(c(140, 160), c(300, 400))
chrs <- chromatogram(msd, rt = rtr, mz = mzr)

## Each row of the returned MChromatograms object corresponds to one mz-rt
## range. The Chromatogram for the first range in the first file is empty,
## because the retention time range is outside of the file's rt range:
chrs[1, 1]
#> Object of class: Chromatogram
#> Intensity values aggregated using: sum
#> length of object: 0
#> from file: 1
#> mz range: [NA, NA]
#> MS level: 1

## The mz and/or rt ranges used are provided as featureData of the object
fData(chrs)
#>   mzmin mzmax rtmin rtmax polarity
#> 1   140   160    10    60        1
#> 2   300   320   280   300        1

## The mz method can be used to extract the m/z ranges directly
mz(chrs)
#>      mzmin mzmax
#> [1,]   140   160
#> [2,]   300   320

## Get the extracted chromatogram for the first range in the second file
chr <- chrs[2, 1]
chr
#> Object of class: Chromatogram
#> Intensity values aggregated using: sum
#> length of object: 148
#> from file: 2
#> mz range: [140.0022, 159.9989]
#> rt range: [10.24602, 59.71602]
#> MS level: 1

plot(rtime(chr), intensity(chr), xlab = "rtime", ylab = "intensity")
