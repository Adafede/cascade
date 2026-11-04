library(tinytest)
library(cascade)

fct_reorder2_base <- utils::getFromNamespace("fct_reorder2_base", "cascade")
dense_rank_base <- utils::getFromNamespace("dense_rank_base", "cascade")
prepare_hierarchy <- utils::getFromNamespace("prepare_hierarchy", "cascade")
hierarchies_progress <- utils::getFromNamespace(
  "hierarchies_progress",
  "cascade"
)
hierarchies_grouped_progress <- utils::getFromNamespace(
  "hierarchies_grouped_progress",
  "cascade"
)
prepare_plot <- utils::getFromNamespace("prepare_plot", "cascade")
prepare_plot_2 <- utils::getFromNamespace("prepare_plot_2", "cascade")
plot_histograms_taxo <- utils::getFromNamespace(
  "plot_histograms_taxo",
  "cascade"
)
plot_histograms_litt <- utils::getFromNamespace(
  "plot_histograms_litt",
  "cascade"
)
treemaps_progress <- utils::getFromNamespace("treemaps_progress", "cascade")
treemaps_progress_no_title <- utils::getFromNamespace(
  "treemaps_progress_no_title",
  "cascade"
)

# --- 1. fct_reorder2_base & dense_rank_base ---
f_unused <- factor(c("a", "b", "a", "b"), levels = c("a", "b", "unused"))
r_desc <- fct_reorder2_base(
  f_unused,
  c(1, 1, 2, 2),
  c(10, 50, 20, 100),
  desc = TRUE
)
expect_identical(levels(r_desc), c("b", "a", "unused"))

r_asc <- fct_reorder2_base(
  f_unused,
  c(1, 1, 2, 2),
  c(10, 50, 20, 100),
  desc = FALSE
)
expect_identical(levels(r_asc), c("a", "b", "unused"))

expect_identical(
  dense_rank_base(c(30, 10, 20, 10, NA_real_, 30)),
  c(3L, 1L, 2L, 1L, NA_integer_, 3L)
)
expect_identical(dense_rank_base(numeric(0L)), integer(0L))

# --- 2. prepare_hierarchy validation, rescale = TRUE vs FALSE, safety == "Y" ---
expect_error(
  prepare_hierarchy(data.frame(), type = "invalid"),
  pattern = "'type' must be either 'analysis' or 'literature'"
)

anal_df <- data.table::data.table(
  id = rep("s1", 6L),
  sample = rep("s1", 6L),
  species = rep("Swertia chirayita", 6L),
  peak_id = 1:6,
  feature_id = 1:6,
  smiles_2D = c("C1", "C2", "C3", NA_character_, "C5", "C6"),
  inchikey_2D = c("IK1", "IK2", "IK3", NA_character_, "IK5", "IK6"),
  score_biological = rep(0.8, 6L),
  score_chemical = rep(0.8, 6L),
  score_final = c(0.9, 0.9, 0.9, NA_real_, 0.2, 0.9),
  best_candidate_organism = rep("Swertia", 6L),
  best_candidate_1 = c(
    "Terpenoids",
    "Alkaloids",
    "Shikimates",
    "notAnnotated",
    "notConfident",
    "notClassified"
  ),
  best_candidate_2 = c(
    "Monoterpenoids",
    "Indole",
    "Flavonoids",
    "notAnnotated notAnnotated",
    "notConfident notConfident",
    "notClassified notClassified"
  ),
  best_candidate_3 = c(
    "Monoterpenoids", # same as best_candidate_2 -> triggers safety == "Y" ('Monoterpenoids sub')
    "Simple indole",
    "Flavones",
    "notAnnotated notAnnotated notAnnotated",
    "notConfident notConfident notConfident",
    "notClassified notClassified notClassified"
  ),
  comparison_score = rep(0.9, 6L),
  intensity = c(100, 200, 300, 50, 50, 50),
  integral = c(10, 20, 30, 5, 5, 5)
)

h_unscaled <- prepare_hierarchy(
  anal_df,
  type = "analysis",
  detector = "ms",
  rescale = FALSE
)
h_rescaled <- prepare_hierarchy(
  anal_df,
  type = "analysis",
  detector = "ms",
  rescale = TRUE
)
expect_true(any(grepl("Monoterpenoids sub", h_unscaled$labels, fixed = TRUE)))
# When rescale = TRUE, sum of values / 3 equals 1
expect_equal(sum(h_rescaled$values) / 3, 1.0, tolerance = 1e-12)
expect_true(sum(h_unscaled$values) > 10)

# --- 3. hierarchies_progress & hierarchies_grouped_progress ---
lit_tbl <- data.table::data.table(
  structure = c("Q1", "Q2", "Q3"),
  chemical_pathway = c("Terpenoids", "Alkaloids", "Shikimates"),
  chemical_superclass = c("Monoterpenoids", "Indole", "Flavonoids"),
  chemical_class = c("Secoiridoids", "Simple", "Flavones"),
  taxaLabels = c("Swertia chirayita", "Swertia perennis", "Swertia var. alba"),
  taxa = c("QT1", "QT2", "QT3"),
  references = c("QR1", "QR2", "QR3"),
  referencesLabels = c("R1", "R2", "R3")
)
lit_single_sp <- data.table::data.table(
  structure = "Q1",
  chemical_pathway = "Terpenoids",
  chemical_superclass = "Monoterpenoids",
  chemical_class = "Secoiridoids",
  taxaLabels = "Swertia chirayita",
  taxa = "QT1",
  references = "QR1",
  referencesLabels = "R1"
)
lit_empty <- data.table::data.table()

# Multi-word comparison keeps full species name in organism/sample
hp_multi <- hierarchies_progress(
  list(Swertia = lit_tbl, Empty = lit_empty),
  comparison = c("Swertia chirayita", "Swertia perennis")
)
expect_true(any(hp_multi$Swertia$species == "Swertia chirayita"))
expect_identical(nrow(hp_multi$Empty), 1L)
expect_identical(hp_multi$Empty$values, 0)
expect_true(is.na(hp_multi$Empty$parents))

# Grouped progress: >1 non-var taxa returns hierarchy; ==1 non-var taxon returns NULL; empty returns 1-row NA
hgp <- hierarchies_grouped_progress(
  list(Multi = lit_tbl, Single = lit_single_sp, Empty = lit_empty)
)
expect_inherits(hgp$Multi, "data.table")
expect_true(nrow(hgp$Multi) > 0L)
expect_null(hgp$Single)
expect_identical(nrow(hgp$Empty), 1L)
expect_identical(hgp$Empty$values, 0)

# --- 4. prepare_plot (empty presamples, without 'Other', and >5 subgroups color overflow fallback) ---
pp_empty <- prepare_plot(
  data.table::data.table(
    parents = "",
    ids = "Terpenoids",
    labels = "Terpenoids",
    values = 10,
    sample = "s1",
    species = "Swertia"
  )
)
expect_identical(nrow(pp_empty), 0L)

# 7 subgroups under 1 parent (exceeds 5 shades in microshades palette -> triggers #00%04d fallback)
overflow_h <- data.table::data.table(
  parents = rep("Terpenoids", 7L),
  ids = paste0("Terpenoids-Sub", 1:7),
  labels = paste0("Sub", 1:7),
  values = seq(70, 10, by = -10),
  sample = rep("s1", 7L),
  species = rep("Swertia", 7L)
)
pp_over <- prepare_plot(overflow_h)
expect_identical(nrow(pp_over), 7L)
expect_true(any(grepl("^#00000", as.character(pp_over$color))))

# Plot with xlab = FALSE
p_litt_no_xlab <- plot_histograms_litt(pp_over, label = "Sample", xlab = FALSE)
expect_inherits(p_litt_no_xlab, "ggplot")

# --- 5. prepare_plot_2 (empty table and all 5 biological score tiers) & plot_histograms_taxo(mode = "neg") ---
pp2_empty <- prepare_plot_2(data.table::data.table(inchikey_2D = NA_character_))
expect_identical(nrow(pp2_empty), 0L)

pp2_tiers_in <- data.table::data.table(
  sample = rep("s1", 5L),
  species = rep("Swertia chirayita", 5L),
  peak_id = 1:5,
  peak_rt_apex = as.numeric(1:5),
  peak_area = rep(100, 5L),
  feature_id = 1:5,
  feature_area = rep(1000, 5L),
  rt = as.numeric(1:5),
  mz = seq(100, 500, by = 100),
  comparison_score = rep(0.9, 5L),
  smiles_2D = paste0("C", 1:5),
  inchikey_2D = paste0("IK", 1:5),
  score_biological = c(0.95, 0.85, 0.65, 0.30, 0.10),
  score_chemical = rep(0.8, 5L),
  score_final = rep(0.8, 5L),
  best_candidate_organism = rep("Swertia", 5L),
  best_candidate_1 = c(
    "Terpenoids",
    "Alkaloids",
    NA_character_,
    "notConfident",
    "notClassified"
  ),
  best_candidate_2 = c(
    "Monoterpenoids",
    NA_character_,
    "notConfident notConfident",
    "Sub4",
    "Sub5"
  ),
  best_candidate_3 = c(
    "Secoiridoids",
    NA_character_,
    "notConfident notConfident notConfident",
    "notClassified",
    "Cla5"
  )
)
pp2_tiers <- prepare_plot_2(pp2_tiers_in)
expect_identical(nrow(pp2_tiers), 5L)
expect_true(all(
  c("Species", "Genus", "Family", "Kingdom", "Other") %in%
    as.character(pp2_tiers$name_2)
))

chrom_trace <- data.table::data.table(
  rtime = seq(0, 6, by = 0.5),
  intensity = seq(1, 13, by = 1)
)
p_taxo_neg <- plot_histograms_taxo(
  pp2_tiers,
  chromatogram = chrom_trace,
  level = "max",
  mode = "neg",
  time_min = 0,
  time_max = 6
)
expect_inherits(p_taxo_neg, "ggplot")

# --- 6. treemaps_progress (regular + "special" comparative domain, with/without structure column) ---
h_swertia <- hp_multi$Swertia
h_gentiana <- data.table::copy(h_swertia)
h_gentiana[, `:=`(sample = "Gentiana", species = "Gentiana")]
h_special <- data.table::rbindlist(list(
  data.table::copy(h_swertia)[, `:=`(sample = "Swertia", species = "Swertia")],
  h_gentiana
))
h_map <- list(Swertia = h_swertia, Gentiana = h_gentiana, special = h_special)
t_map <- list(
  Swertia = data.table::data.table(structure = c("Q1", "Q1", "Q2")),
  Gentiana = data.table::data.table(other_col = 1:4)
)

tm_titled <- treemaps_progress(
  xs = c("Swertia", "Gentiana", "special"),
  type = "treemap",
  hierarchies = h_map,
  tables = t_map
)
expect_identical(names(tm_titled), c("Swertia", "Gentiana", "special"))
expect_inherits(tm_titled$Swertia, "plotly")
expect_inherits(tm_titled$special, "plotly")

tm_null_tables <- treemaps_progress(
  xs = "Swertia",
  type = "sunburst",
  hierarchies = h_map,
  tables = NULL
)
expect_inherits(tm_null_tables$Swertia, "plotly")
