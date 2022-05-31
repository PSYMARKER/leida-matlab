function Plot_K_tpm(data_dir,selectedK)
%
% Plot the mean transition probability matrices for each condition.
%
% INPUT:
% data_dir      directory with the results from running the hypothesis
%               tests on the state-to-state transition probabilities
% selectedK     K defined by the user
%
% OUTPUT:
% .fig/.png     plot of the estimated mean TPM for each condition
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with results for the transition probabilities (output from LEiDA_stats_TransitionMatrix.m)
file_TM = 'LEiDA_Stats_TransitionMatrix.mat';

% Load required data:
load([data_dir file_TM],'cond','TMnorm','Index_Conditions');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Code below generates the figures relative to the transition probabilities
% Code below generates the figures relative to the fractional occupancy
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITION PROBABILITIES FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Mean transition matrix for each conditionPLOT OF THE MEAN TRANSITION MATRIX FOR EACH CONDITION
TM_cond = cell(1,n_Cond);
for j = 1:n_Cond
    % Compute mean TPM omitting the NaN values
    TM_cond{j} = squeeze(mean(TMnorm(Index_Conditions == j,:,:),'omitnan'));
end

% Plot of mean TPM
disp(' ');
disp(['Plotting the mean transition probability matrix with ' num2str(selectedK) ' PL states for each condition:'])
for i = 1:n_Cond
    Fig = figure('Position', get(0, 'Screensize'));
    mean_TM = TM_cond{i};
    imagesc(mean_TM, 'AlphaData', 0.8);
    for c_out = 1:selectedK
        for c_in = 1:selectedK
            caption = sprintf('%.2f', mean_TM(c_out,c_in));
            text(c_in,c_out, caption, 'FontSize', 15, 'Color', [1, 1, 1],'FontWeight','bold','HorizontalAlignment','Center');
        end
    end
    colormap('jet');
    set(colorbar,'YTick',[0 0.2 0.4 0.6 0.8 1], 'ylim', [0 1], 'FontSize', 12);
    ticklabels = cell(1,selectedK);
    for p = 1:selectedK
        ticklabels{p} = p;
    end
    xticks(1:1:selectedK); xticklabels(ticklabels);
    yticks(1:1:selectedK); yticklabels(ticklabels);
    set(gca,'FontSize',14,'FontWeight','bold')
    axis square
    ylabel('From FC State')
    xlabel('To FC State')
    
    saveas(Fig, fullfile(data_dir, ['K' num2str(selectedK) '_meanTPM_' cond{i} '.png']),'png');
    saveas(Fig, fullfile(data_dir, ['K' num2str(selectedK) '_meanTPM_' cond{i} '.fig']),'fig');
    disp(['- Mean TPM of condition ' cond{i} ' successfully saved as K' num2str(selectedK) '_meanTPM_' cond{i}]);
end

close all;