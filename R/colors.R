library(microshades, quietly = TRUE)

## Colors
nice_colors <- list(
  microshades_palette("micro_cvd_green", lightest = FALSE),
  microshades_palette("micro_cvd_orange", lightest = FALSE),
  microshades_palette("micro_cvd_blue", lightest = FALSE),
  microshades_palette("micro_cvd_purple", lightest = FALSE),
  # microshades_palette("micro_cvd_gray", lightest = FALSE),
  microshades_palette("micro_cvd_turquoise", lightest = FALSE),
  microshades_palette("micro_purple", lightest = FALSE),
  microshades_palette("micro_orange", lightest = FALSE),
  microshades_palette("micro_blue", lightest = FALSE),
  microshades_palette("micro_green", lightest = FALSE),
  microshades_palette("micro_brown", lightest = FALSE)
  # microshades_palette("micro_gray", lightest = FALSE)
)

## Grey
grey_colors <- list(
  microshades_palette("micro_cvd_gray", lightest = FALSE),
  microshades_palette("micro_gray", lightest = FALSE)
)

sunburst_colors <- character()

for (i in seq_len(length(nice_colors))) {
  sunburst_colors[i] <- nice_colors[[i]][5]
}

sunburst_colors_grey <- character()

for (i in seq_len(length(grey_colors))) {
  sunburst_colors_grey[i] <- grey_colors[[i]][5]
}

## From microshades
cvd_green <- c("#4E7705", "#6D9F06", "#97CE2F", "#BDEC6F", "#DDFFA0")
cvd_orange <- c("#9D654C", "#C17754", "#F09163", "#FCB076", "#FFD5AF")
cvd_blue <- c("#098BD9", "#56B4E9", "#7DCCFF", "#BCE1FF", "#E7F4FF")
cvd_turquoise <- c("#148F77", "#009E73", "#43BA8F", "#48C9B0", "#A3E4D7")
cvd_purple <- c("#7D3560", "#A1527F", "#CC79A7", "#E794C1", "#EFB6D6")
orange <- c("#ff7f00", "#fe9929", "#fdae6b", "#fec44f", "#feeda0")
purple <- c("#6a51a3", "#807dba", "#9e9ac8", "#bcbddc", "#dadaeb")

## Qualitative
paired <- c(
  "#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#ffff99",
  "#b15928"
)

minibits <- c(
  "#817",
  "#a35",
  "#c66",
  "#e94",
  "#ed0",
  "#9d5",
  "#4d8",
  "#2cb",
  "#0bc",
  "#09c",
  "#36b",
  "#639"
)

tableau20 <- c(
  "#4E79A7",
  "#A0CBE8",
  "#F28E2B",
  "#FFBE7D",
  "#59A14F",
  "#8CD17D",
  "#B6992D",
  "#F1CE63",
  "#499894",
  "#86BCB6",
  "#E15759",
  "#FF9D9A",
  "#79706E",
  "#BAB0AC",
  "#D37295",
  "#FABFD2",
  "#B07AA1",
  "#D4A6C8",
  "#9D7660",
  "#D7B5A6"
)

elife <- c(
  "#7CB13F",
  "#336A2D",
  "#2994D2",
  "#08589B",
  "#D71D62",
  "#861450"
)

#' Sequential
green_4 <- c(
  "#e5efd9",
  "#b0d08c",
  "#7cb13f",
  "#4a6a26"
)

blue_4 <- c(
  "#d4eaf6",
  "#7fbfe4",
  "#2994d2",
  "#19597e"
)

pink_4 <- c(
  "#f7d2e0",
  "#e777a1",
  "#d71d62",
  "#81113b"
)

green_24 <- c(
  "#fafcf7",
  "#eff6e8",
  "#e5efd9",
  "#dae9c9",
  "#d0e3ba",
  "#c5ddab",
  "#bbd69b",
  "#b0d08c",
  "#a6ca7c",
  "#9bc46d",
  "#91bd5e",
  "#86b74e",
  "#7cb13f",
  "#72a33a",
  "#689535",
  "#5e8730",
  "#54782b",
  "#4a6a26",
  "#405c21",
  "#374e1c",
  "#2d4017",
  "#233212",
  "#19230d",
  "#0f1508"
)

blue_24 <- c(
  "#f6fbfd",
  "#e5f2fa",
  "#d4eaf6",
  "#c3e1f2",
  "#b2d8ef",
  "#a1d0eb",
  "#90c7e8",
  "#7fbfe4",
  "#6db6e0",
  "#5caedd",
  "#4ba5d9",
  "#3a9dd6",
  "#2994d2",
  "#2688c1",
  "#227cb0",
  "#1f70a0",
  "#1c658f",
  "#19597e",
  "#154d6d",
  "#12415c",
  "#0f354c",
  "#0b293b",
  "#081e2a",
  "#051219"
)

pink_24 <- c(
  "#fdf6f9",
  "#fae4ec",
  "#f7d2e0",
  "#f4c0d3",
  "#f1aec6",
  "#ed9cba",
  "#ea89ad",
  "#e777a1",
  "#e46594",
  "#e15388",
  "#dd417b",
  "#da2f6f",
  "#d71d62",
  "#c61b5a",
  "#b51852",
  "#a3164a",
  "#921443",
  "#81113b",
  "#700f33",
  "#5f0d2b",
  "#4d0a23",
  "#3c081b",
  "#2b0614",
  "#1a030c"
)
