# Development and assessment of analytic methods to improve the measurement of cognition in longitudinal studies of aging through the use of sub-studies with comprehensive neuropsychological testing
Code accompaniment to the paper, the goal of the analysis was to use data from the Harmonized Assessment Protocol Survey to both develop new approaches to estimate summary scores of cognitive functioning and evaluate the performance of these approaches both against each other and against existing approaches. 
We evaluate four approaches (two novel, two existing): 
1. CFA (novel): we use fixed parameter calibration to estimate measurement models in the ELSA-HCAP sample and apply them to ELSA
2. Regression-based prediction (novel): We use regression models estimated in the ELSA-HCAP sample to predict cognition in the full ELSA study
3. Word recall (existing): Use word recall as a measure of cognitive functioning
4. Summary score (existing): Mean of z-scores from assessments of memory, orientation, and verbal fluency
The available code files provide code for core data analysis steps in the paper:
- 01-estimateCFA.R: Estimate CFA models to create gold standard cognitive factor scores and the full CFA for prediction to the ELSA sample
- 02-prediction.R: Estimate regression model to apply to the full ELSA sample
- 03-applytoELSA.R: Calculate novel measures in the full ELSA sample
- 04-bootstrap-ELSAHCAP.R: Calculate bootstrapped estimates of performance in ELSA-HCAP using different train/test sets

