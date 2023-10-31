# crankTurning-dynamicPrimatives

This repository contains the code used in the writing of the paper: 

James Hermus, Joeseph Doeringer, Dagmar Sternad, and Neville Hogan, "Dynamic Primitives in Constrained Action: Systematic Changes in the Zero-Force Trajectory" Currently accepted to the Journal of Neurophysiology. [Link](https://journals.physiology.org/doi/abs/10.1152/jn.00082.2023) 

## Code Organization and Data
The experimental data file is quite large and not stored on Git Hub. Download the [crank_data](https://www.dropbox.com/sh/goa68d93k81vd98/AADflYlnum4B5ZK0uSr05dssa?dl=0) folder from dropbox and place it in the crankTurning-dynamicPrimitives folder. The ReadMe.txt file within the crank_data folder provides further details about the data set itself.

[main](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/main.m): This script generates the figures presented in the paper. Processing all trials, conditions, and stiffness levels takes some time. For this reason, first a check is performed to determine if the analysis.mat data structure already exists. If so it will load it. If not it will generate the analysis.mat.

The code has two main classes:

- [trialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/trialAnalysis.m): This class imports and processes each trial
- [crossTrialAnalysis](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/crossTrialAnalysis.m): This class combines information from multiple trials and conditions. This class also contains the statistical analysis and plotting funtions. 

Lastly, there is [lnCoordinatesExampleFig](https://github.com/jameshermus/crankTurning-dynamicPrimitives/blob/main/lnCoordinatesExampleFig.m): This function generates the conceptual Figure 9 for the paper.

