INPUT=$1
sirius \
--input $1 \
--output ${INPUT%.*} \
config \
--IsotopeSettings.filter true \
--FormulaSearchDB BIO,COCONUT \
--Timeout.secondsPerTree 10 \
--FormulaSettings.enforced HCNOPS \
--Timeout.secondsPerInstance 10 \
--AdductSettings.detectable "[M]-, [M + Cl]-, [M + Br]-, [M + CH2O2 - H]-, [M - H2 + Na]-, [M + C2H4O2 - H]-, [M + C2HF3O2 - H]-, [M + C2H3N - H]-, [M - 2H + K ]-, [M - H2O - H]-, [M - H]-, [M - H4O2 - H]-, [M - H6O3 - H]-, [M - H8O4 - H]-, [M - H10O5 - H]-, [M - CO - H]-, [M - C2O2 - H]-, [M - C3O3 - H]-, [M - C2H4 - H]-, [M - C2H5 - H]-, [M - CH2O - H]-, [M - CH4O - H]-, [M - C2H2O - H]-, [M - CO2 - H]-, [M - C2O4 - H]-, [M - CH2O2 - H]-, [M - C2H6O - H]-, [M - C2H4O2 - H]-, [M - CH6O3 - H]-, [M - C2H2O3 - H]-, [M - C3H7O2 - H]-, [M - C3H4O3 - H]-, [M - C4H8O2 - H]-, [M - H2O4S - H]-, [M - H3O4P - H]-, [M - C5H10O2 - H]-, [M - C3H4O4 - H]-, [M - C5H8O4 - H]-, [M - C9H6O2 - H]-, [M - C6H10O4 - H]-, [M - C12H20O8 - H]-, [M - C9H8O2 - H]-, [M - C7H4O4 - H]-, [M - C9H6O3 - H]-, [M - C6H10O5 - H]-, [M - C12H20O10 - H]-, [M - C18H30O15 - H]-" \
--UseHeuristic.mzToUseHeuristicOnly 650 \
--AlgorithmProfile orbitrap \
--IsotopeMs2Settings IGNORE \
--MS2MassDeviation.allowedMassDeviation 5.0ppm \
--NumberOfCandidatesPerIon 1 \
--UseHeuristic.mzToUseHeuristic 300 \
--FormulaSettings.detectable B,Br,Cl,F,I,Se,Si \
--NumberOfCandidates 10 \
--ZodiacNumberOfConsideredCandidatesAt300Mz 10 \
--ZodiacRunInTwoSteps true \
--ZodiacEdgeFilterThresholds.minLocalConnections 10 \
--ZodiacEdgeFilterThresholds.thresholdFilter 0.95 \
--ZodiacEpochs.burnInPeriod 2000 \
--ZodiacEpochs.numberOfMarkovChains 10 \
--ZodiacNumberOfConsideredCandidatesAt800Mz 50 \
--ZodiacEpochs.iterations 20000 \
--AdductSettings.enforced , \
--AdductSettings.fallback "[M]-, [M - H2 + Na]-, [M + C2H4O2 - H]-, [M + C2HF3O2 - H]-, [M - H2O - H]-, [M - H]-" \
--FormulaResultThreshold true \
--InjectElGordoCompounds true \
--StructureSearchDB BIO,COCONUT \
--RecomputeResults false \
formula \
zodiac \
fingerprint \
structure \
canopus \
W