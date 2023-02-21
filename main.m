% Data analysis for "Dynamic Primitives in Constrained Action: Systematic 
% Changes in the Zero-Force Trajectory"
% Filename:	main.m
% Author:  James Hermus
% Date:    21 Feb 2023
% Description:	This script generates the figures presented in the paper

clear all
close all
clc

dbstop if error
restoredefaultpath;

subjDex = 1:10; % number of subjects
speedDex = 1:6; % speed and direction number
trialDex = 1:21; % trial number

%% Check crank_data is downloaded
if(7~=exist('crank_data', 'dir'))
    error('The crank_data folder is missing. Make sure to download crank_data folder from dropbox and add it to your local repository. Information is provided in the README.md.');
end

%% Load already created analysis structures
tic
filee = ['analysis.mat'];
if(exist(filee)>1) % If file already exists load
    load(filee);
else
    parfor stiffDex = 1:9
        analysis{stiffDex} = crossTrialAnalysis(subjDex,speedDex,trialDex,stiffDex);
    end
    save('-v7.3','analysis');
end
toc

%% Figure 3: This figure was replicated from Hermus et al. 2020, "Separating
%% neural influences from peripheral mechanics: the speed-curvature relation
%% in mechanically constrained actions"

%% Figure 4: ZFT representative subject single trial
subj_rep = 2;
trial_rep = 16;
analysis{5}.test{subj_rep,1,trial_rep}.get_zftPlot();
analysis{5}.test{subj_rep,2,trial_rep}.get_zftPlot();
analysis{5}.test{subj_rep,3,trial_rep}.get_zftPlot();

analysis{5}.test{subj_rep,4,trial_rep}.get_zftPlot();
analysis{5}.test{subj_rep,5,trial_rep+3}.get_zftPlot();
analysis{5}.test{subj_rep,6,trial_rep}.get_zftPlot();

%% Figure 5: ZFT all subjects
analysis{5}.get_subjectAveZFTplots();

%% Figure 6: ln(r) PC plots
anovaPlot = true;
subjectPlot = false;
analysis{5}.get_lnrCoordinate_plots(subjectPlot,anovaPlot);

% Check other stiffness levels and generate analysis.outputt_lnr_PC
subjectPlot = false;
anovaPlot = false;
for stiffDex = [1:4,6:9]    
    analysis{stiffDex}.get_lnrCoordinate_plots(subjectPlot,anovaPlot);
    close
end

%% Figure 7: ln(r) PC plots as stiffness varies
figure; hold on;
for stiffDex = 1:9
    analysis{stiffDex}.get_lnrCoordinates_plots_allStiff(); hold on;
end

%% Figure 8: CV of V and statistics
analysis{5}.get_CVplots();
analysis{5}.outputt_cv_v.P

%% Figure 9: Insight to ln(r) coordinates
lnCoordinatesExampleFig();

