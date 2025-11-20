

<!-- README.md is generated from README.qmd. Please edit that file -->

# CASCADE <img src='https://raw.githubusercontent.com/adafede/cascade/main/man/figures/logo.svg' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![No Maintenance
Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)
[![CRAN
status](https://www.r-pkg.org/badges/version/tima.png)](https://CRAN.R-project.org/package=cascade)
[![R-CMD-check](https://github.com/adafede/cascade/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adafede/cascade/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/adafede/cascade/graph/badge.svg)](https://app.codecov.io/gh/adafede/cascade)
[![r-universe
badge](https://adafede.r-universe.dev/cascade/badges/version?&color=blue&style=classic.png)](https://adafede.r-universe.dev/cascade)
<!-- [![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/adafede/cascade/) -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14515158.svg)](https://doi.org/10.5281/zenodo.14515158)
<!-- badges: end -->

**C**ontextualizing untargeted **A**nnotation with **S**emi-quantitative
**C**harged **A**erosol **D**etection for pertinent characterization of
natural **E**xtracts.

⚠️ This repository is not maintained and will not be except for extreme
interest. It has just been opened for the sake of transparency.

The initial work is available at
https://doi.org/10.1021/acs.jafc.3c03099, with some improvements made
since then. The workflow is illustrated below.

![Workflow](https://raw.githubusercontent.com/adafede/cascade/main/man/figures/cascade-workflow.svg)  

## Requirements

Here is what you *minimally* need:

- A file
  ([.mzML](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#mzML))
  containing [DDA
  MS](https://en.wikipedia.org/w/index.php?title=Data-dependent_acquisition)
  data with an additional detector
  ([PDA](https://en.wikipedia.org/wiki/UV_detectors),
  [ELSD](https://en.wikipedia.org/wiki/Evaporative_light_scattering_detector),
  [CAD](https://en.wikipedia.org/wiki/Charged_aerosol_detector))
  - In case you don’t know how to obtain it, see:
    [wiki/How-to-create-a-compliant-mzML-file](https://github.com/Adafede/cascade/wiki/How-to-create-a-compliant-mzML-file)
- A file ([.csv](https://en.wikipedia.org/wiki/Comma-separated_values))
  containing features, as obtained by [mzmine](https://mzio.io)
- A file ([.tsv](https://en.wikipedia.org/wiki/Tab-separated_values))
  containing annotations, as obtained by
  [TIMA](https://taxonomicallyinformedannotation.github.io/tima)

## Installation

As the package is not (yet) available on CRAN, you will need to install
with:

``` r
install.packages(
  "cascade",
  repos = c(
    "https://adafede.r-universe.dev",
    "https://bioc.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
```

Once installed, you are ready to go through our
[documentation](https://adafede.github.io/cascade/articles), with the
major steps detailed.

## Main Citations

### CASCADE

<https://doi.org/10.1021/acs.jafc.3c03099>

According to which steps you used, please give credit to the authors of
the tools/resources used.

### mzmine

General: <https://doi.org/10.1038/s41587-023-01690-2>

### SIRIUS

General: <https://doi.org/10.1038/s41592-019-0344-8>

- *CSI:FingerId*: <https://doi.org/10.1073/pnas.1509788112>
- *ZODIAC*: <https://doi.org/10.1038/s42256-020-00234-6>
- *CANOPUS*: <https://doi.org/10.1038/s41587-020-0740-8>
- *COSMIC*: <https://doi.org/10.1038/s41587-021-01045-9>

### LOTUS

General: <https://doi.org/10.7554/eLife.70780>

⚠️ Do not forget to cite which version you used:
<https://doi.org/10.5281/zenodo.5794106>

### ISDB

General: <https://doi.org/10.1021/acs.analchem.5b04804>

⚠️ Do not forget to cite which version you used:
<https://doi.org/10.5281/zenodo.5607185>

### TIMA

General: <https://doi.org/10.3389/fpls.2019.01329>

⚠️ Do not forget to cite which version you used:
<https://doi.org/10.5281/zenodo.5797920>

### Others

- NPClassifier: <https://doi.org/10.1021/acs.jnatprod.1c00399>

## Additional software credits

| Package | Version | Citation |
|:---|:---|:---|
| base | 4.5.2 | R Core Team (2025) |
| baseline | 1.3.7 | Liland, Almøy, and Mevik (2010); Liland and Mevik (2025) |
| BiocManager | 1.30.27 | Morgan and Ramos (2025) |
| BiocParallel | 1.44.0 | Wang et al. (2025) |
| BiocVersion | 3.22.0 | Morgan (2025) |
| cascade | 0.0.0.9001 | Rutz and Wolfender (2023); Rutz (2025) |
| caTools | 1.18.3 | Tuszynski (2024) |
| data.table | 1.17.8 | Barrett et al. (2025) |
| gt | 1.1.0 | Iannone et al. (2025) |
| htmltools | 0.5.8.1 | Cheng et al. (2024) |
| knitr | 1.50 | Xie (2014); Xie (2015); Xie (2025) |
| MSnbase | 2.36.0 | Gatto and Lilley (2012); Gatto, Gibb, and Rainer (2020) |
| mzR | 2.44.0 | Pedrioli et al. (2004); Keller et al. (2005); Kessner et al. (2008); Martens et al. (2010); Chambers et al. (2012) |
| pkgload | 1.4.1 | Wickham et al. (2025) |
| plotly | 4.11.0 | Sievert (2020) |
| R.utils | 2.13.0 | Bengtsson (2025) |
| rmarkdown | 2.30 | Xie, Allaire, and Grolemund (2018); Xie, Dervieux, and Riederer (2020); Allaire et al. (2025) |
| stringi | 1.8.7 | Gagolewski (2022) |
| testthat | 3.3.0 | Wickham (2011) |
| tidytable | 0.11.2 | Fairbanks (2024) |
| tidyverse | 2.0.0 | Wickham et al. (2019) |
| tima | 2.12.0 | Rutz et al. (2019); Rutz and Allard (2025) |
| WikidataQueryServiceR | 1.0.0 | Popov (2020) |

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-rmarkdown2025" class="csl-entry">

Allaire, JJ, Yihui Xie, Christophe Dervieux, Jonathan McPherson, Javier
Luraschi, Kevin Ushey, Aron Atkins, et al. 2025.
*<span class="nocase">rmarkdown</span>: Dynamic Documents for r*.
<https://github.com/rstudio/rmarkdown>.

</div>

<div id="ref-datatable" class="csl-entry">

Barrett, Tyson, Matt Dowle, Arun Srinivasan, Jan Gorecki, Michael
Chirico, Toby Hocking, Benjamin Schwendinger, and Ivan Krylov. 2025.
*<span class="nocase">data.table</span>: Extension of
“<span class="nocase">data.frame</span>”*. <https://r-datatable.com>.

</div>

<div id="ref-Rutils" class="csl-entry">

Bengtsson, Henrik. 2025. *<span class="nocase">R.utils</span>: Various
Programming Utilities*. <https://henrikbengtsson.github.io/R.utils/>.

</div>

<div id="ref-mzR2012" class="csl-entry">

Chambers, Matthew C., Maclean, Brendan, Burke, Robert, Amodei, et al.
2012. “<span class="nocase">A cross-platform toolkit for mass
spectrometry and proteomics</span>.” *Nat Biotech* 30 (10): 918–20.
<https://doi.org/10.1038/nbt.2377>.

</div>

<div id="ref-htmltools" class="csl-entry">

Cheng, Joe, Carson Sievert, Barret Schloerke, Winston Chang, Yihui Xie,
and Jeff Allen. 2024. *<span class="nocase">htmltools</span>: Tools for
HTML*. <https://github.com/rstudio/htmltools>.

</div>

<div id="ref-tidytable" class="csl-entry">

Fairbanks, Mark. 2024. *<span class="nocase">tidytable</span>: Tidy
Interface to “<span class="nocase">data.table</span>”*.
<https://markfairbanks.github.io/tidytable/>.

</div>

<div id="ref-stringi" class="csl-entry">

Gagolewski, Marek. 2022. “<span class="nocase">stringi</span>: Fast and
Portable Character String Processing in R.” *Journal of Statistical
Software* 103 (2): 1–59. <https://doi.org/10.18637/jss.v103.i02>.

</div>

<div id="ref-MSnbase2020" class="csl-entry">

Gatto, Laurent, Sebastian Gibb, and Johannes Rainer. 2020. “MSnbase,
Efficient and Elegant r-Based Processing and Visualisation of Raw Mass
Spectrometry Data.” *bioRxiv*.

</div>

<div id="ref-MSnbase2012" class="csl-entry">

Gatto, Laurent, and Kathryn Lilley. 2012. “MSnbase - an r/Bioconductor
Package for Isobaric Tagged Mass Spectrometry Data Visualization,
Processing and Quantitation.” *Bioinformatics* 28: 288–89.

</div>

<div id="ref-gt" class="csl-entry">

Iannone, Richard, Joe Cheng, Barret Schloerke, Ellis Hughes, Alexandra
Lauer, JooYoung Seo, Ken Brevoort, and Olivier Roy. 2025.
*<span class="nocase">gt</span>: Easily Create Presentation-Ready
Display Tables*. <https://gt.rstudio.com>.

</div>

<div id="ref-mzR2005" class="csl-entry">

Keller, Andrew, Jimmy Eng, Ning Zhang, Xiao-jun Li, and Ruedi Aebersold.
2005. “A Uniform Proteomics MS/MS Analysis Platform Utilizing Open XML
File Formats.” *Mol Syst Biol*.

</div>

<div id="ref-mzR2008" class="csl-entry">

Kessner, Darren, Matt Chambers, Robert Burke, David Agus, and Parag
Mallick. 2008. “ProteoWizard: Open Source Software for Rapid Proteomics
Tools Development.” *Bioinformatics* 24 (21): 2534–36.
<https://doi.org/10.1093/bioinformatics/btn323>.

</div>

<div id="ref-baseline2010" class="csl-entry">

Liland, Kristian Hovde, Trygve Almøy, and Bjørn-Helge Mevik. 2010.
“Optimal Choice of Baseline Correction for Multivariate Calibration of
Spectra.” *Applied Spectroscopy* 64: 1007–16.
<https://doi.org/10.1366/000370210792434350>.

</div>

<div id="ref-baseline2025" class="csl-entry">

Liland, Kristian Hovde, and Bjørn-Helge Mevik. 2025.
*<span class="nocase">baseline</span>: Baseline Correction of Spectra*.
<https://github.com/khliland/baseline/>.

</div>

<div id="ref-mzR2010" class="csl-entry">

Martens, Lennart, Matthew Chambers, Marc Sturm, Darren Kessner, Fredrik
Levander, Jim Shofstahl, Wilfred H Tang, et al. 2010. “MzML - a
Community Standard for Mass Spectrometry Data.” *Mol Cell Proteomics*.
<https://doi.org/10.1074/mcp.R110.000133>.

</div>

<div id="ref-BiocVersion" class="csl-entry">

Morgan, Martin. 2025. *BiocVersion: Set the Appropriate Version of
Bioconductor Packages*. <https://doi.org/10.18129/B9.bioc.BiocVersion>.

</div>

<div id="ref-BiocManager" class="csl-entry">

Morgan, Martin, and Marcel Ramos. 2025. *BiocManager: Access the
Bioconductor Project Package Repository*.
<https://doi.org/10.32614/CRAN.package.BiocManager>.

</div>

<div id="ref-mzR2004" class="csl-entry">

Pedrioli, Patrick G A, Jimmy K Eng, Robert Hubley, Mathijs Vogelzang,
Eric W Deutsch, Brian Raught, Brian Pratt, et al. 2004. “A Common Open
Representation of Mass Spectrometry Data and Its Application to
Proteomics Research.” *Nat Biotechnol* 22 (11): 1459–66.
<https://doi.org/10.1038/nbt1031>.

</div>

<div id="ref-WikidataQueryServiceR" class="csl-entry">

Popov, Mikhail. 2020. *WikidataQueryServiceR: API Client Library for
“Wikidata Query Service”*.
<https://github.com/bearloga/WikidataQueryServiceR>.

</div>

<div id="ref-base" class="csl-entry">

R Core Team. 2025. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

<div id="ref-cascade2025" class="csl-entry">

Rutz, Adriano. 2025. *<span class="nocase">cascade</span>:
Contextualizing Untargeted Annotation with Semi-Quantitative Charged
Aerosol Detection for Pertinent Characterization of Natural Extracts*.

</div>

<div id="ref-tima2025" class="csl-entry">

Rutz, Adriano, and Pierre-Marie Allard. 2025.
*<span class="nocase">tima</span>: Taxonomically Informed Metabolite
Annotation*. <https://doi.org/10.5281/zenodo.5797920>.

</div>

<div id="ref-tima2019" class="csl-entry">

Rutz, Adriano, Miwa Dounoue-Kubo, Simon Ollivier, Jonathan Bisson,
Mohsen Bagheri, Tongchai Saesong, Samad Nejad Ebrahimi, Kornkanok
Ingkaninan, Jean-Luc Wolfender, and Pierre-Marie Allard. 2019.
“Taxonomically Informed Scoring Enhances Confidence in Natural Products
Annotation.” *Frontiers in Plant Science* 10.
<https://doi.org/10.3389/FPLS.2019.01329>.

</div>

<div id="ref-cascade2023" class="csl-entry">

Rutz, Adriano, and Jean-Luc Wolfender. 2023. “Automated Composition
Assessment of Natural Extracts: Untargeted Mass Spectrometry-Based
Metabolite Profiling Integrating Semiquantitative Detection.” *Journal
of Agricultural and Food Chemistry* 71 (46).
<https://doi.org/10.1021/acs.jafc.3c03099>.

</div>

<div id="ref-plotly" class="csl-entry">

Sievert, Carson. 2020. *Interactive Web-Based Data Visualization with r,
Plotly, and Shiny*. Chapman; Hall/CRC. <https://plotly-r.com>.

</div>

<div id="ref-caTools" class="csl-entry">

Tuszynski, Jarek. 2024. *<span class="nocase">caTools</span>: Tools:
Moving Window Statistics, GIF, Base64, ROC AUC, Etc*.

</div>

<div id="ref-BiocParallel" class="csl-entry">

Wang, Jiefei, Martin Morgan, Valerie Obenchain, Michel Lang, Ryan
Thompson, and Nitesh Turaga. 2025. *BiocParallel: Bioconductor
Facilities for Parallel Evaluation*.
<https://doi.org/10.18129/B9.bioc.BiocParallel>.

</div>

<div id="ref-testthat" class="csl-entry">

Wickham, Hadley. 2011. “<span class="nocase">testthat</span>: Get
Started with Testing.” *The R Journal* 3: 5–10.
<https://journal.r-project.org/articles/RJ-2011-002/>.

</div>

<div id="ref-tidyverse" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

</div>

<div id="ref-pkgload" class="csl-entry">

Wickham, Hadley, Winston Chang, Jim Hester, and Lionel Henry. 2025.
*<span class="nocase">pkgload</span>: Simulate Package Installation and
Attach*. <https://github.com/r-lib/pkgload>.

</div>

<div id="ref-knitr2014" class="csl-entry">

Xie, Yihui. 2014. “<span class="nocase">knitr</span>: A Comprehensive
Tool for Reproducible Research in R.” In *Implementing Reproducible
Computational Research*, edited by Victoria Stodden, Friedrich Leisch,
and Roger D. Peng. Chapman; Hall/CRC.

</div>

<div id="ref-knitr2015" class="csl-entry">

———. 2015. *Dynamic Documents with R and Knitr*. 2nd ed. Boca Raton,
Florida: Chapman; Hall/CRC. <https://yihui.org/knitr/>.

</div>

<div id="ref-knitr2025" class="csl-entry">

———. 2025. *<span class="nocase">knitr</span>: A General-Purpose Package
for Dynamic Report Generation in R*. <https://yihui.org/knitr/>.

</div>

<div id="ref-rmarkdown2018" class="csl-entry">

Xie, Yihui, J. J. Allaire, and Garrett Grolemund. 2018. *R Markdown: The
Definitive Guide*. Boca Raton, Florida: Chapman; Hall/CRC.
<https://bookdown.org/yihui/rmarkdown>.

</div>

<div id="ref-rmarkdown2020" class="csl-entry">

Xie, Yihui, Christophe Dervieux, and Emily Riederer. 2020. *R Markdown
Cookbook*. Boca Raton, Florida: Chapman; Hall/CRC.
<https://bookdown.org/yihui/rmarkdown-cookbook>.

</div>

</div>
