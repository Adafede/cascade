library(tinytest)
library(cascade)

y_as_na <- utils::getFromNamespace("y_as_na", "cascade")
dt_left_join <- utils::getFromNamespace("dt_left_join", "cascade")
dt_inner_join <- utils::getFromNamespace("dt_inner_join", "cascade")
dt_full_join <- utils::getFromNamespace("dt_full_join", "cascade")
p_acn_i <- utils::getFromNamespace("p_acn_i", "cascade")
predict_response <- utils::getFromNamespace("predict_response", "cascade")
molinfo <- utils::getFromNamespace("molinfo", "cascade")
keep_best_candidates <- utils::getFromNamespace(
  "keep_best_candidates",
  "cascade"
)
make_confident <- utils::getFromNamespace("make_confident", "cascade")
make_other <- utils::getFromNamespace("make_other", "cascade")
no_other <- utils::getFromNamespace("no_other", "cascade")
make_no_stereo <- utils::getFromNamespace("make_no_stereo", "cascade")
make_chromatographiable <- utils::getFromNamespace(
  "make_chromatographiable",
  "cascade"
)
prepare_comparison <- utils::getFromNamespace("prepare_comparison", "cascade")

# --- 1. y_as_na (factor and character inputs) ---
expect_identical(
  y_as_na(factor(c("a", "", "b")), ""),
  c("a", NA_character_, "b")
)
expect_identical(
  y_as_na(c("x", "none", NA_character_), "none"),
  c("x", NA_character_, NA_character_)
)

# --- 2. dt_left_join, dt_inner_join, dt_full_join (named by, colliding cols,
# suffixes) ---
dx <- data.table::data.table(
  id_x = c("b", "a", "c"),
  shared = c(1L, 1L, 1L),
  val = c(10, 20, 30)
)
dy <- data.table::data.table(
  id_y = c("a", "b", "d"),
  shared = c(1L, 1L, 1L),
  val = c(100, 200, 400),
  extra = c("A", "B", "D")
)

# Named by with mixed named/unnamed elements + colliding column 'val'
lj <- dt_left_join(
  dx,
  dy,
  by = c("id_x" = "id_y", "shared"),
  suffix = c("_left", "_right")
)
expect_identical(
  names(lj),
  c("id_x", "shared", "val_left", "val_right", "extra")
)
expect_identical(lj$id_x, c("b", "a", "c"))
expect_identical(lj$val_left, c(10, 20, 30))
expect_identical(lj$val_right, c(200, 100, NA_real_))
expect_identical(lj$extra, c("B", "A", NA_character_))

ij <- dt_inner_join(
  dx,
  dy,
  by = c("id_x" = "id_y", "shared"),
  suffix = c("_left", "_right")
)
expect_identical(ij$id_x, c("a", "b"))
expect_identical(ij$val_left, c(20, 10))
expect_identical(ij$val_right, c(100, 200))

fj <- dt_full_join(
  dx,
  dy,
  by = c("id_x" = "id_y", "shared"),
  suffix = c("_left", "_right")
)
expect_identical(nrow(fj), 4L)
expect_identical(
  names(fj),
  c("id_x", "shared", "val_left", "val_right", "extra")
)
expect_identical(fj$id_x, c("b", "a", "c", "d"))

# Natural join (by = NULL)
nat_lj <- dt_left_join(
  data.table::data.table(k = c(2L, 1L), x = c("b", "a")),
  data.table::data.table(k = c(1L, 2L), y = c("A", "B")),
  by = NULL
)
expect_identical(nat_lj$k, c(2L, 1L))
expect_identical(nat_lj$y, c("B", "A"))

# --- 3. p_acn_i & predict_response & molinfo ---
expect_equal(
  p_acn_i(c(0, 1, 2), q1 = 2, q2 = -3, q3 = 5),
  c(5, 4, 7),
  tolerance = 1e-12
)
expect_identical(
  predict_response(acn = c(20, 50, 80), peak_area = c(-5, 0, NA_real_)),
  c(0, 0, 0)
)
# Hand-verified calculation at acn = 50, peak_area = 100
p1 <- 0.00001 * 2500 - 0.0006 * 50 - 0.0778
p2 <- 0.00002 * 2500 - 0.00022 * 50 + 0.05499
p3 <- -0.00017 * 2500 + 0.0209 * 50 + 1.4041
den <- 2 * p1 + p2
expected_resp <- 100^(1 / den) * exp(-p3 / den)
expect_equal(
  predict_response(acn = 50, peak_area = 100),
  expected_resp,
  tolerance = 1e-12
)

mi <- molinfo("CCO")
expect_true(grepl(
  "https://www.simolecule.com/cdkdepict/depict/bot/svg\\?smi=CCO",
  mi
))

# --- 4. keep_best_candidates (with and without rank_final, pipes, sentinels)
# ---
raw_cand <- data.table::data.table(
  feature_id = c(1L, 1L, 2L, 3L, 4L),
  feature_mz = c("300.1|300.1", "300.1", "400.2", "500.3", "600.4"),
  feature_rt = c("5.0", "5.0", "6.0", "7.0", "8.0"),
  candidate_structure_smiles_no_stereo = c("C=C|C=O", "C=O", "", "CCC", "CCCC"),
  candidate_structure_inchikey_connectivity_layer = c(
    "IK1|IK2",
    "IK2",
    "",
    "IK3",
    "IK4"
  ),
  candidate_structure_molecular_formula = c(
    "C2H4",
    "CH2O",
    "",
    "C3H8",
    "C4H10"
  ),
  candidate_count_similarity_peaks_matched = c("5", "4", "", "3", "6"),
  candidate_score_similarity = c("0.9", "0.8", "", "0.7", "0.95"),
  score_initial = c("0.9", "0.8", "", "0.7", "0.95"),
  score_biological = c("0.5", "0.9", "", "0.5", "0.8"),
  score_chemical = c("0.5", "0.9", "", "0.5", "0.8"),
  score_final = c("0.4", "0.95", "", "0.6", "0.85"),
  candidate_structure_organism_occurrence_closest = c(
    "Swertia",
    "Swertia",
    "",
    "Gentiana",
    "Swertia"
  ),
  feature_pred_tax_npc_01pat_val = c(
    "Terpenoids",
    "Terpenoids",
    "Alkaloids",
    "Shikimates",
    "Polyketides"
  ),
  feature_pred_tax_npc_02sup_val = c(
    "Mono",
    "Mono",
    "Indole",
    "Flavonoids",
    "Macrolides"
  ),
  feature_pred_tax_npc_03cla_val = c(
    "Iridoids",
    "Iridoids",
    "Simple",
    "Flavones",
    "Ansamycins"
  ),
  feature_pred_tax_npc_01pat_score = c("0.9", "0.9", "0.8", "0.7", "0.9"),
  feature_pred_tax_npc_02sup_score = c("0.9", "0.9", "0.8", "0.7", "0.9"),
  feature_pred_tax_npc_03cla_score = c("0.9", "0.9", "0.8", "0.7", "0.9"),
  candidate_structure_tax_npc_01pat = c(
    "Terpenoids or Polyketides",
    "Alkaloids or Other",
    "",
    "",
    "Terpenoids"
  ),
  candidate_structure_tax_npc_02sup = c(
    "Monoterpenoids",
    "Indole or Other",
    "",
    "",
    ""
  ),
  candidate_structure_tax_npc_03cla = c(
    "Secoiridoids",
    "Simple indole",
    "",
    "",
    ""
  )
)
# Without rank_final column: feature 1 should pick row 2 (score_final = 0.95 >
# 0.4)
kbc_no_rk <- keep_best_candidates(raw_cand)
expect_identical(nrow(kbc_no_rk), 4L)
expect_identical(kbc_no_rk[feature_id == "1"]$smiles_2D, "C=O")
expect_identical(kbc_no_rk[feature_id == "1"]$best_candidate_1, "Alkaloids")
expect_identical(kbc_no_rk[feature_id == "1"]$best_candidate_2, "Indole")
# Feature 2 has empty smiles -> notAnnotated
expect_identical(kbc_no_rk[feature_id == "2"]$best_candidate_1, "notAnnotated")
expect_identical(
  kbc_no_rk[feature_id == "2"]$best_candidate_2,
  "notAnnotated notAnnotated"
)
expect_identical(
  kbc_no_rk[feature_id == "2"]$best_candidate_3,
  "notAnnotated notAnnotated notAnnotated"
)
# Feature 3 has smiles but empty classes -> notClassified
expect_identical(kbc_no_rk[feature_id == "3"]$best_candidate_1, "notClassified")
expect_identical(
  kbc_no_rk[feature_id == "3"]$best_candidate_2,
  "notClassified notClassified"
)
expect_identical(
  kbc_no_rk[feature_id == "3"]$best_candidate_3,
  "notClassified notClassified notClassified"
)
# Feature 4 has best_candidate_1 = Terpenoids, but empty 2 & 3
expect_identical(kbc_no_rk[feature_id == "4"]$best_candidate_1, "Terpenoids")
expect_identical(
  kbc_no_rk[feature_id == "4"]$best_candidate_2,
  "Terpenoids notClassified"
)
expect_identical(
  kbc_no_rk[feature_id == "4"]$best_candidate_3,
  "Terpenoids notClassified notClassified"
)

# With rank_final column: rank_final = 1 overrides lower score_final
raw_cand_rk <- data.table::copy(raw_cand)
raw_cand_rk[, rank_final := c("1", "2", "1", "1", "1")]
kbc_rk <- keep_best_candidates(raw_cand_rk)
expect_identical(kbc_rk[feature_id == "1"]$smiles_2D, "C=C")
expect_identical(kbc_rk[feature_id == "1"]$best_candidate_1, "Terpenoids")

# --- 5. make_confident exact boundary ---
mc_in <- data.table::data.table(
  score_final = c("0.5", "0.4999", NA_character_),
  best_candidate_1 = c("A", "B", "C"),
  best_candidate_2 = c("A2", "B2", "C2"),
  best_candidate_3 = c("A3", "B3", "C3")
)
mc_out <- make_confident(mc_in, score = 0.5)
expect_identical(mc_out$best_candidate_1, c("A", "notConfident", NA_character_))
expect_identical(
  mc_out$best_candidate_2,
  c("A2", "notConfident notConfident", NA_character_)
)
expect_identical(
  mc_out$best_candidate_3,
  c("A3", "notConfident notConfident notConfident", NA_character_)
)

# --- 6. make_other & no_other ---
# Empty table branch
empty_mo <- make_other(data.table::data.table())
expect_identical(nrow(empty_mo), 0L)

# Top-4 superclass retention per pathway, 5th+ collapsed into "Other"
mo_in <- data.table::data.table(
  smiles_2D = paste0("S", 1:6),
  inchikey_2D = paste0("IK", 1:6),
  best_candidate_1 = rep("Terpenoids", 6L),
  best_candidate_2 = c("Sup1", "Sup2", "Sup3", "Sup4", "Sup5", "Sup6"),
  best_candidate_3 = paste0("Cla", 1:6),
  custom_area = c(100, 90, 80, 70, 20, 10)
)
mo_out <- make_other(mo_in, value = "custom_area")
expect_identical(
  mo_out[order(-custom_area)]$best_candidate_2,
  c("Sup1", "Sup2", "Sup3", "Sup4", "Other", "Other")
)

no_in <- data.table::data.table(
  best_candidate_1 = c(
    "Alkaloids",
    "notAnnotated",
    "notClassified",
    "notConfident",
    NA_character_
  )
)
expect_identical(no_other(no_in)$best_candidate_1, "Alkaloids")

# --- 7. make_no_stereo & make_chromatographiable ---
ns_in <- data.table::data.table(
  taxon = c("T1", "T1", "T1"),
  structure = c("Q1", "Q2", "Q3"),
  structure_smiles_2D = c("C=C", "C=C", "C=O"),
  reference_doi = c("10.1/a", "10.1/a", "10.1/a")
)
ns_out <- make_no_stereo(ns_in)
expect_identical(nrow(ns_out), 2L)
expect_false("structure_smiles_2D" %in% names(ns_out))

mc_struct <- data.table::data.table(
  id = 1:9,
  structure_exact_mass = c(50, 49.999, 1500, 1500.001, 500, 500, 500, 500, 500),
  structure_xlogp = c(2, 2, 2, 2, -1, -1.001, 6, 6.001, 0)
)
mc_filt <- make_chromatographiable(
  mc_struct,
  mass_min = 50,
  mass_max = 1500,
  logp_min = -1,
  logp_max = 6
)
expect_identical(mc_filt$id, c(1L, 3L, 5L, 7L, 9L))

# --- 8. prepare_comparison (taxonomic prioritization, thresholds, ACN
# correction, empty precor) ---
local({
  tmp_inf <- tempfile(fileext = ".tsv")
  tmp_ninf <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_inf, tmp_ninf)), add = TRUE)

  # Peak 1: feature 1 has exact species match (Swertia chirayita), feature 2
  # has genus match (Swertia) -> only feature 1 kept
  # Peak 2: feature 3 has genus match (Swertia), feature 4 has unrelated
  # organism (Gentiana) -> only feature 3 kept
  # Peak 3: feature 5 & 6 have unrelated organisms -> both kept (sum == 0)
  # Peak 4: feature 7 has comparison_score = 0.5 < min_similarity_prefilter
  # (0.6) -> moved to peaks_min_precor
  inf_dt <- data.table::data.table(
    sample = rep("s1", 7L),
    peak_id = c(1L, 1L, 2L, 2L, 3L, 3L, 4L),
    peak_rt_min = c(1, 1, 2, 2, 3, 3, 4),
    peak_rt_apex = c(1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5),
    peak_rt_max = c(2, 2, 3, 3, 4, 4, 5),
    peak_area = rep(100, 7L),
    feature_id = 1:7,
    feature_rt = c(1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5),
    feature_mz = seq(100, 700, by = 100),
    feature_area = rep(1000, 7L),
    comparison_score = c(0.9, 0.85, 0.6, 0.9, 0.8, 0.88, 0.5)
  )
  ninf_dt <- data.table::data.table(
    sample = "s1",
    peak_id = NA_integer_,
    peak_rt_min = NA_real_,
    peak_rt_apex = NA_real_,
    peak_rt_max = NA_real_,
    peak_area = NA_real_,
    feature_id = 8L,
    feature_rt = 9.0,
    feature_mz = 800,
    feature_area = 500,
    comparison_score = NA_real_
  )
  data.table::fwrite(inf_dt, tmp_inf, sep = "\t")
  data.table::fwrite(ninf_dt, tmp_ninf, sep = "\t")

  cand_conf <- data.table::data.table(
    feature_id = 1:8,
    rt = c(1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5, 9.0),
    mz = seq(100, 800, by = 100),
    smiles_2D = paste0("C", 1:8),
    inchikey_2D = paste0("IK", 1:8),
    best_candidate_1 = rep("Terpenoids", 8L),
    best_candidate_2 = rep("Monoterpenoids", 8L),
    best_candidate_3 = rep("Secoiridoids", 8L),
    best_candidate_organism = c(
      "Swertia chirayita",
      "Swertia",
      "Swertia",
      "Gentiana",
      "Gentiana",
      "Arabidopsis",
      "Swertia chirayita",
      "Swertia chirayita"
    ),
    species = rep("Swertia chirayita", 8L),
    score_biological = rep(0.9, 8L),
    score_chemical = rep(0.9, 8L),
    score_final = rep(0.9, 8L),
    mode = rep("pos", 8L)
  )

  comp_res <- prepare_comparison(
    features_informed = tmp_inf,
    features_not_informed = tmp_ninf,
    candidates_confident = cand_conf,
    min_similarity_prefilter = 0.6,
    min_similarity_filter = 0.8,
    mode = "pos",
    correct_acn = TRUE
  )
  expect_identical(comp_res$peaks_maj_precor$feature_id, as.numeric(1:6))
  # After taxonomic prioritization: feature 1 (peak 1), feature 3 (peak 2),
  # features 5 & 6 (peak 3)
  expect_identical(comp_res$peaks_maj_precor_taxo$feature_id, c(1, 3, 5, 6))
  # After min_similarity_filter >= 0.8: feature 3 (score 0.75) is filtered out
  expect_identical(comp_res$peaks_maj_precor_taxo_cor$feature_id, c(1, 5, 6))
  # ACN correction modified peak_area from 100
  expect_true(all(comp_res$peaks_maj$peak_area != 100))

  # Empty peaks_maj_precor branch (min_similarity_prefilter = 0.99)
  comp_empty <- prepare_comparison(
    features_informed = tmp_inf,
    features_not_informed = tmp_ninf,
    candidates_confident = cand_conf,
    min_similarity_prefilter = 0.99,
    min_similarity_filter = 0.99,
    mode = "pos",
    correct_acn = FALSE
  )
  expect_identical(nrow(comp_empty$peaks_maj_precor_taxo), 0L)
  expect_identical(nrow(comp_empty$peaks_maj_precor_taxo_cor), 0L)
})
