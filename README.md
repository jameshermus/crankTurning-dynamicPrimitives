# crankTurning-dynamicPrimatives

This repository contains the code use in the writing of the paper: 

James Hermus, Joeseph Doeringer, Dagmar Sternad, and Neville Hogan, "Dynamic Primitives in Constrained Action: Systematic Changes in the Zero-Force Trajectory" Currently submitted to the Joural of Neurophysiology. 

## Code Organization and Data
The experimental data file is quite large and not stored on github. Download the [crank_data](https://www.dropbox.com/sh/goa68d93k81vd98/AADflYlnum4B5ZK0uSr05dssa?dl=0) folder from dropbox and place it in the crankTurning-dynamicPrimitives folder. The ReadMe.txt file with in the crank_data folder provids futher details about the data set its self.

[main](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/main.m): This script generates the figures presented in the paper. Processing all trials, condistions, and stiffness levels takes some time. For this reason, first a check is preformed to determine if the analysis.mat data structure already exists. If so it will load it. If not it will generate the analysis.mat.

The code has two main classes:

- [trialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/trialAnalysis.m): This class imports and processes each trial
- [crossTrialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/crossTrialAnalysis.m): This class combines information from multiple trials and conditions. This class also contains the statisticall analysis and plotting funtions. 

Lastly, there is [lnCoordinatesExampleFig](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/lnCoordinatesExampleFig.m): This function generates the conceptual Figure 9 for the paper.

