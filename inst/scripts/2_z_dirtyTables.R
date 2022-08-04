cat(
  "Features:",
  nrow(annotations |>
    dplyr::filter(mode == "pos")),
  "\n",
  "cadPeaks:",
  nrow(df_meta_cad_pos |>
    dplyr::filter(peak_id != 0)),
  "\n",
  "FeaturesLinked:",
  nrow(compared_peaks_list_cad_pos$peaks_all),
  "\n",
  "featuresPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(featuresPerPeak_old != 0) |>
      dplyr::pull(featuresPerPeak_old)
  ),
  "\n",
  "structuresPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(structuresPerPeak_old != 0) |>
      dplyr::pull(structuresPerPeak_old)
  ),
  "\n",
  "MFPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(molecularFormulasPerPeak_old != 0) |>
      dplyr::pull(molecularFormulasPerPeak_old)
  ),
  "\n",
  "classesPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(chemicalClassesPerPeak_old != 0) |>
      dplyr::pull(chemicalClassesPerPeak_old)
  ),
  "\n"
)

cat(
  "Features:",
  nrow(annotations |>
    dplyr::filter(mode == "neg")),
  "\n",
  "cadPeaks:",
  nrow(df_meta_cad_neg |>
    dplyr::filter(peak_id != 0)),
  "\n",
  "FeaturesLinked:",
  nrow(compared_peaks_list_cad_neg$peaks_all),
  "\n",
  "featuresPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(featuresPerPeak_old != 0) |>
      dplyr::pull(featuresPerPeak_old)
  ),
  "\n",
  "structuresPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(structuresPerPeak_old != 0) |>
      dplyr::pull(structuresPerPeak_old)
  ),
  "\n",
  "MFPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(molecularFormulasPerPeak_old != 0) |>
      dplyr::pull(molecularFormulasPerPeak_old)
  ),
  "\n",
  "classesPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(chemicalClassesPerPeak_old != 0) |>
      dplyr::pull(chemicalClassesPerPeak_old)
  ),
  "\n"
)

cat(
  "correlatedfeatures",
  nrow(
    compared_peaks_list_cad_pos$peaks_maj_precor |>
      dplyr::filter(comparison_score >= params$chromato$peak$similarity$filter)
  ),
  "\n",
  "correlatedfeaturesPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedFeaturesPerPeak_old != 0) |>
      dplyr::pull(correlatedFeaturesPerPeak_old)
  ),
  "\n",
  "correlatedstructuresPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedStructuresPerPeak_old != 0) |>
      dplyr::pull(correlatedStructuresPerPeak_old)
  ),
  "\n",
  "correlatedMFPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedMolecularFormulasPerPeak_old != 0) |>
      dplyr::pull(correlatedMolecularFormulasPerPeak_old)
  ),
  "\n",
  "correlatedclassesPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedChemicalClassesPerPeak_old != 0) |>
      dplyr::pull(correlatedChemicalClassesPerPeak_old)
  ),
  "\n"
)

cat(
  "correlatedFeatures",
  nrow(
    compared_peaks_list_cad_neg$peaks_maj_precor |>
      dplyr::filter(comparison_score >= params$chromato$peak$similarity$filter)
  ),
  "\n",
  "correlatedfeaturesPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedFeaturesPerPeak_old != 0) |>
      dplyr::pull(correlatedFeaturesPerPeak_old)
  ),
  "\n",
  "correlatedstructuresPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedStructuresPerPeak_old != 0) |>
      dplyr::pull(correlatedStructuresPerPeak_old)
  ),
  "\n",
  "correlatedMFPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedMolecularFormulasPerPeak_old != 0) |>
      dplyr::pull(correlatedMolecularFormulasPerPeak_old)
  ),
  "\n",
  "correlatedclassesPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(correlatedChemicalClassesPerPeak_old != 0) |>
      dplyr::pull(correlatedChemicalClassesPerPeak_old)
  ),
  "\n"
)

cat(
  "finalfeatures",
  nrow(
    compared_peaks_list_cad_pos$peaks_maj_precor_taxo_cor |>
      dplyr::filter(comparison_score >= params$chromato$peak$similarity$filter)
  ),
  "\n",
  "finalfeaturesPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalFeaturesPerPeak_old != 0) |>
      dplyr::pull(finalFeaturesPerPeak_old)
  ),
  "\n",
  "finalstructuresPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalStructuresPerPeak_old != 0) |>
      dplyr::pull(finalStructuresPerPeak_old)
  ),
  "\n",
  "finalMFPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalMolecularFormulasPerPeak_old != 0) |>
      dplyr::pull(finalMolecularFormulasPerPeak_old)
  ),
  "\n",
  "finalclassesPerPeak:",
  mean(
    df_meta_cad_pos |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalChemicalClassesPerPeak_old != 0) |>
      dplyr::pull(finalChemicalClassesPerPeak_old)
  ),
  "\n"
)

cat(
  "finalfeatures",
  nrow(
    compared_peaks_list_cad_neg$peaks_maj_precor_taxo_cor |>
      dplyr::filter(comparison_score >= params$chromato$peak$similarity$filter)
  ),
  "\n",
  "finalfeaturesPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalFeaturesPerPeak_old != 0) |>
      dplyr::pull(finalFeaturesPerPeak_old)
  ),
  "\n",
  "finalstructuresPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalStructuresPerPeak_old != 0) |>
      dplyr::pull(finalStructuresPerPeak_old)
  ),
  "\n",
  "finalMFPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalMolecularFormulasPerPeak_old != 0) |>
      dplyr::pull(finalMolecularFormulasPerPeak_old)
  ),
  "\n",
  "finalclassesPerPeak:",
  mean(
    df_meta_cad_neg |>
      dplyr::filter(peak_id != 0) |>
      dplyr::filter(finalChemicalClassesPerPeak_old != 0) |>
      dplyr::pull(finalChemicalClassesPerPeak_old)
  ),
  "\n"
)
