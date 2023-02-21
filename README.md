# crankTurning-dynamicPrimatives

This repository contains the code use in the writting of the paper: 

James Hermus, Joeseph Doeringer, Dagmar Sternad, and Neville Hogan, "Dynamic Primitives in Constrained Action: Systematic Changes in the Zero-Force Trajectory" Currently submitted to the Joural of Neurophysiology. 

## Code Organization and Data
The experimenal data file is quite large and not stored on github. Download the [crank_data](https://www.dropbox.com/scl/fo/tiybowr0f0vkru7o85gij/h?dl=0&rlkey=udwtmb0zrnlm8atrk1w6sc2uc) folder from dropbox and place it in the crankTurning-dynamicPrimitives folder.

The code is organized into three sections:

- [trialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/trialAnalysis.m): This class imports and processes each trial
- [crossTrialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/crossTrialAnalysis.m): This class combines information from multiple trials and conditions. This class also contains the statisticall analysis and plotting funtions. 
- [lnCoordinatesExampleFig](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/lnCoordinatesExampleFig.m): This class generates Figure 9 from the paper
- [main](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/main.m): This script generates the figures presented in the paper. Processing all trials, condistions, and stiffness levels takes some time. For this reason, first a check is preformed to determine if the analysis.mat data structure already exists. If so it will load it. If not it will generate the analysis.mat.
