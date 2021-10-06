;
  wdt:P225 ?taxon_name.
  ?structure wdt:P235 ?structure_id;
  wdt:P233 ?structureSmiles;
  p:P703 ?statement.
  ?statement ps:P703 ?taxon;
  prov:wasDerivedFrom ?ref.
  ?ref pr:P248 ?art.
  ?art wdt:P1476 ?art_title;
  wdt:P356 ?art_doi;
  wdt:P577 ?art_date.
  FILTER(((YEAR(?art_date)) >= 