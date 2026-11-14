# Getting help with cascade

## Where to check first

1. Read the package overview in [README.md](../README.md) and the vignettes at
   <https://adafede.github.io/cascade/articles/>.
2. Check open and closed issues at <https://github.com/adafede/cascade/issues>
   to see whether the problem has already been reported or fixed.

## Filing an issue

If you find a bug or need help with a workflow failure, open an issue at
<https://github.com/adafede/cascade/issues/new> and include:

- The exact function call and arguments that failed.
- Whether the problem occurs with the bundled example data
  (`show_example = TRUE`) or custom `.mzML`, `.csv`, and `.tsv` inputs.
- For custom inputs, which detectors are present in the `.mzML` file (`CAD`,
  `PDA`, `ELSD`, MS1) and the full error message or warning output.
- Your R and package versions (`sessionInfo()`).
