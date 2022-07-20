# CASCADE

**C**ontextualizing *in-depth* untargeted **A**nnotation with **S**emi-quantitative **C**harged **A**erosol **D**etection for pertinent characterization of natural **E**xtracts.

## Get example file (TODO)

## Run peakonly

```
python python/peakonly.py -i inst/extdata/source/mzml/ -o inst/extdata/interim/peakonly/210705_AR_06_V_03_2_01_peakonly.tsv
```

## Run peakonly2mzmine

```
Rscript inst/scripts/peakonly2mzmine.R
```

## Run MZmine

```
bash bash/create_batch_file.sh inst/extdata/source/mzml/qcmix/191113_AR_QCmix_07_09_Pos.mzML
```

```
mzmine --batch config/params/batch_191113_AR_QCmix_07_09_Pos_mzmine.xml
```