sirius \
--output $1 \
config \
--IsotopeSettings.filter true \
--FormulaSearchDB BIO,COCONUT \
--Timeout.secondsPerTree 0 \
--FormulaSettings.enforced HCNOPS \
--Timeout.secondsPerInstance 0 \
--AdductSettings.detectable "[M - H4O2 + H]+, [M + H]+, [M + H3N + H]+, [M + Na]+, [M + K]+, [M - H2O + H]+" \
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
--AdductSettings.fallback "[M + H]+, [M + Na]+, [M + K]+" \
--FormulaResultThreshold true \
--InjectElGordoCompounds true \
--StructureSearchDB BIO,COCONUT \
--RecomputeResults false \
formula \
zodiac \
fingerprint \
structure \
canopus