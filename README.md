# CASCADE

<!-- badges: start -->
[![Lifecycle:Maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
<!-- badges: end -->

**C**ontextualizing untargeted **A**nnotation with **S**emi-quantitative **C**harged **A**erosol **D**etection for pertinent characterization of natural **E**xtracts.

:warning: This repository is not functional, not maintained and will not be except for extreme interest.
 
It has just been opened for the sake of transparency.


1. Run `bash bash/create_batch_file.sh <PATH_TO_YOUR_FILE> <polarity>`

2. Run `mzmine -b <PATH_TO_BATCH_CREATED_IN_1>`

3. Run `Rscript inst/scripts/1_process_compare_peaks.R`

4. Run annotation (TODO)

5. Run `Rscript inst/scripts/2_process_plot_pseudochromatgrams.R`

...

