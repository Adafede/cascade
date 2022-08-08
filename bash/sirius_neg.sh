sirius \
--output $1 \
config \
--IsotopeSettings.filter true \
--FormulaSearchDB BIO,COCONUT \
--Timeout.secondsPerTree 0 \
--FormulaSettings.enforced HCNOPS \
--Timeout.secondsPerInstance 0 \
--AdductSettings.detectable "[M + Cl]-, [M + Br]-, [M - H2O - H]-, [M - H]-" \
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
--AdductSettings.fallback "[M]-, [M + Cl]-, [M + Br]-, [M + CH2O2 - H]-, [M - H + Na - H]-, [M + C2H4O2 - H]-, [M + C2HF3O2 - H]-, [M + C2H3N - H]-, [M - H + K - H]-, [M - H2O - H]-, [M - H]-" \
--FormulaResultThreshold true \
--InjectElGordoCompounds true \
--StructureSearchDB BIO,COCONUT \
--RecomputeResults false \
formula \
zodiac \
fingerprint \
structure \
canopus