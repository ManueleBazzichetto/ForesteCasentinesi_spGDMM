# spGDMM Foreste Casentinesi

This repo includes the R scripts to run [spGDMM models](https://doi.org/10.1111/2041-210X.14259) using vegetation data from the Foreste Casentinesi.

DataFromEVA.R: extract vegetation data for the Foreste Casentinesi from the EVA dataset;
ClimTopoData.R: get climatic, topographic and disturbance data at the location of the vegetation plots;
DataForGDM: clean-up vegetation data and prepare input objects for the spGDMM.
ComputeClimIndices: function for computing a set of climatic indices using the ClimInd R package.