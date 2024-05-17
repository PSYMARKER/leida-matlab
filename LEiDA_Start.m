function LEiDA_Start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function runs LEiDA on any given dataset.
%
% Since we do not know a priori the optimal number of PL states that
% differentiate between conditions, this function analyses the solutions
% obtained for different numbers of PL states (K).
%
% After analysing the output figures saved in the folder LEiDA_Results, the
% user can then choose the optimal number of PL states for subsequent
% detailed analyses using the functions LEiDA_AnalysisK.m and
% LEiDA_AnalysisCentroid.m.
%
% This function contains two sections:
%       (A) User defines the parameters and properties of the study.
%       (B) Run LEiDA and statistics.
%       (C) Generate and save figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis:
%    - Compute the leading eigenvectors for all participants
%    - Cluster the leading eigenvectors of all participants
%    - Compute statistics to compare across conditions
% C: Figures
%    - Analysis of Fractional Occupancy values
%    - Analysis of Dwell Time values
%    - Pyramid of PL states
%
% Tutorial: README.md
% Version:  V1.0, June 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

%% A: STUDY PARAMETERS

% Directory of the LEiDA toolbox folder:
LEiDA_directory = '/Users/vaniamiguel/Documents/leida-matlab-master';
% Directory of the folder with the parcellated neuroimaging data:
Data_directory = '/Users/vaniamiguel/Documents/leida-matlab-master/Time_series_folder/';
% Name of the run to be used to create the folder to save the data:
run_name = 'Tutorial_Results';
% Tag of conditions given in the parcellated image files:
Conditions_tag = {'CONT','AUT','ASP','PDD_NOS'};
% Parcellation applied to the imaging data (see tutorial):
Parcellation = 'AAL120';
% Number of brain areas to consider for analysis:
N_areas = 94;
% Repetition time (TR) of the fMRI data (if unknown set to 1):
TR = 1;
% Maximum number of TRs for all fMRI sessions:
Tmax = 315;
% Apply temporal filtering to data (0: no; 1: yes)
apply_filter = 0;
% Lowpass frequency of filter (default 0.1):
flp = 0.1;
% Highpass frequency of filter (default 0.01):
fhi = 0.01;

% For the statistics:
% Choose 0 (unpaired) if subjects in different conditions are not the
% same; or 1 (paired) if subjects are the same across conditions.
Paired_tests = 0;
% Number of permutations. For the first analysis to be relatively quick,
% run around 500 permutations, but then increase to 10000 to increase the
% reliability of the final statistical results (p-values) for publication.
n_permutations = 500;
% Number of bootstrap samples within each permutation. For the first
% analysis to be relatively quick, choose around 10, but then increase to
% 500 for more reliable final results.
n_bootstraps = 10;

% For the figure of the pyramid of PL states:
% Direction to plot the FC states/brain ('SideView' or 'TopView'):
CortexDirection = 'SideView';


% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Add the LEiDA_directory to the matlab path
addpath(genpath(LEiDA_directory))

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Create a directory to store the results from the current LEiDA run
if ~exist([LEiDA_directory 'res_' run_name '/'], 'dir')
    mkdir([LEiDA_directory 'res_' run_name '/']);
end
leida_res = [LEiDA_directory 'res_' run_name '/'];

%% B: RUN LEADING EIGENVECTOR DYNAMICS ANALYSIS
% Compute the leading eigenvectors of the data
LEiDA_data(Data_directory,leida_res,N_areas,Tmax,apply_filter,flp,fhi,TR);

% Cluster the leading eigenvectors of all subjects
LEiDA_cluster(leida_res);

% Compute the fractional occupancy and perform hypothesis tests
LEiDA_stats_FracOccup(leida_res,Conditions_tag,Paired_tests,n_permutations,n_bootstraps);

% Compute the dwell time and perform hypothesis tests
LEiDA_stats_DwellTime(leida_res,Conditions_tag,Paired_tests,TR,n_permutations,n_bootstraps);

%% C: MAKE FIGURES

% Generate and save the p-value and barplot plots for fractional occupancy
Plot_FracOccup(leida_res)

% Generate and save the p-value and barplot plots for dwell time
Plot_DwellTime(leida_res)

% Plot the centroids obtained using LEiDA and their overlap with Yeo nets
Plot_Centroid_Pyramid(leida_res,Conditions_tag,Parcellation,N_areas,CortexDirection)