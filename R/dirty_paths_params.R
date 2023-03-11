#' Paths
ANNOTATIONS <- list.files(
  path = file.path(
    paths$data$interim$annotations,
    params$annotation$tool
  ),
  pattern = params$filename$mzml,
  full.names = TRUE,
  recursive = TRUE
)
ANNOTATIONS <- list.files(
  path = file.path(
    paths$data$interim$annotations,
    params$annotation$tool
  ),
  pattern = paste0(
    params$filename$annotation$tima$pos,
    "|",
    params$filename$annotation$tima$neg
  ),
  full.names = TRUE,
  recursive = TRUE
)
EXPORT_DIR <- paths$data$interim$peaks
EXPORT_FILE_BPI <-
  paste(params$filename$mzml,
    "featuresInformed",
    "bpi.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_BPI_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "bpi.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_CAD <-
  paste(params$filename$mzml,
    "featuresInformed",
    "cad.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_CAD_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "cad.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_PDA <-
  paste(params$filename$mzml,
    "featuresInformed",
    "pda.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_PDA_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "pda.tsv.gz",
    sep = "_"
  )
IMPORT_FILE_BPI <-
  list.files(EXPORT_DIR, pattern = "featuresInformed_bpi.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresInformed_bpi.tsv.gz")
  )]
IMPORT_FILE_BPI_2 <-
  list.files(EXPORT_DIR, pattern = "featuresNotInformed_bpi.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresNotInformed_bpi.tsv.gz")
  )]
IMPORT_FILE_CAD <-
  list.files(EXPORT_DIR, pattern = "featuresInformed_cad.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresInformed_cad.tsv.gz")
  )]
IMPORT_FILE_CAD_2 <-
  list.files(EXPORT_DIR, pattern = "featuresNotInformed_cad.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresNotInformed_cad.tsv.gz")
  )]
IMPORT_FILE_PDA <-
  list.files(EXPORT_DIR, pattern = "featuresInformed_pda.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresInformed_pda.tsv.gz")
  )]
IMPORT_FILE_PDA_2 <-
  list.files(EXPORT_DIR, pattern = "featuresNotInformed_pda.tsv.gz")[grepl(
    pattern = gsub(
      pattern = "_.*$",
      replacement = "",
      x = params$filename$mzml
    ),
    x = list.files(EXPORT_DIR, pattern = "featuresNotInformed_pda.tsv.gz")
  )]

# FEATURES <- "~/Downloads/test_quant.csv"
FEATURES <-
  paste0(
    "data/interim/mzmine/lists/",
    params$filename$mzml,
    "_gnps_quant.csv"
  )


# GNPS_JOB <- "97d7c50031a84b9ba2747e883a5907cd"

# TOYSET <- "~/data/20210701_10043/fractions"
# TOYSET <- "~/../../Volumes/LaCie/data/20210701_10043/fractions"
# TOYSET <- "~/data/20210701_10043/test"
TOYSET <- paths$data$source$mzml$path

#' Generic parameters
WORKERS <- params$workers

#' Parameters for LC alignment
TIME_MIN <- params$chromato$time$min
TIME_MAX <- params$chromato$time$max
CAD_SHIFT <- params$chromato$shift$cad
PDA_SHIFT <- params$chromato$shift$pda
ESTIMATED_SOLUBLITIY_LIMIT <- params$misc$solubility$limit

#' Parameters for signal improvement
FOURRIER_COMPONENTS <- params$signal$fourrier$components
FREQUENCY <- params$signal$frequency
RESAMPLE <- params$signal$resample

#' Parameters adapted from Excel sheet from paper shortDOI: 10/ghmvhz
sigma <- params$signal$sigma
k2 <- sigma / params$signal$k2 # 30
k4 <- sigma / params$signal$k4 # 200
smoothing_width <- params$signal$smoothing
baseline_adjust <- params$signal$baseline

#' Parameters related to MS/CAD
INTENSITY_MS_MIN <- params$chromato$intensity$ms1$min
PEAK_SIMILARITY <- params$chromato$peak$similarity$filter
PEAK_SIMILARITY_PREFILTER <-
  params$chromato$peak$similarity$prefilter
RT_TOL <- params$chromato$peak$tolerance$rt
PPM <- params$chromato$peak$tolerance$ppm
AREA_MIN <- params$chromato$peak$area$min

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- params$annotation$confidence$min
