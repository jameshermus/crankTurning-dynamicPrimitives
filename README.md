# crankTurning-dynamicPrimatives

Publication

James Hermus, Joeseph Doeringer, Dagmar Sternad, and Neville Hogan, "Dynamic Primitives in Constrained Action: Systematic Changes in the Zero-Force Trajectory" Currently submitted to the Joural of Neurophysiology. 

Code Organization

The code is organized into three sections:

The experimenal data file is quite large and not stored on github. Download the "crank_data" folder (avalible here) and place it in the crankTurning-dynamicPrimitives folder.

trialAnalysis: This class imports and processes each trial
crossTrialAnalysis: This class combines information from multiple trials and conditions. This class also contains the statisticall analysis and plotting funtions. 
lnCoordinatesExampleFig: This class generates Figure 9 from the paper
main: This script generates the figures presented in the paper. Processing all trials, condistions, and stiffness levels takes some time. For this reason, first a check is preformed to determine if the analysis.mat data structure already exists. If so it will load it. If not it will generate the analysis.mat.


