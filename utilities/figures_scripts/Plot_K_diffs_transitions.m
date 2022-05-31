function Plot_K_diffs_transitions(data_dir,selectedK)
%
% Plot the the results from the hypothesis tests obtained from comparing
% the mean state-to-state transition probability between conditions.

% INPUT:
% data_dir      directory with the results from running the hypothesis
%               tests on the state-to-state transition probabilities
% selectedK     K defined by the user
%
% OUTPUT:
% .fig/.png     plot of the two-sided p-values obtained from the comparison
%               of the mean state-to-state transition probabilities for 
%               each pair of conditions
%
% Authors: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% File with results for the fractional occupancy  (output from LEiDA_stats_FracOccup.m)
file_TM = 'LEiDA_Stats_TransitionMatrix.mat';

% Load required data:
load([data_dir file_TM],'cond','TM_pval2sided');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Create a map to place the figures
pos = power(n_Cond-1,2);
index_fig = reshape(1:pos, n_Cond-1, n_Cond-1).';
subplot_map = ones(n_Cond-1);
subplot_map = triu(subplot_map).';

% Define indices for conditions in subplot
subplot_indices = find(subplot_map);

% Possible pairs of conditions comparisons
condRow = zeros(1,n_Cond*(n_Cond-1)/2);
condCol = zeros(1,n_Cond*(n_Cond-1)/2);
cond_pair = 1;
% disp('Possible pairs of condition comparisons:')
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        % disp([num2str(cond{cond1}) ' - ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

% Plot of differences between transition probabilities across conditions
disp(' ');
disp('Plotting the differences in transition probabilites between conditions:')
Fig = figure('Position', get(0, 'Screensize'));
for s_ind = 1:length(subplot_indices)
    
    % disp(['- ' cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}])
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    % Something in the background
    I = magic(selectedK);
    % Generate where each text will go
    [X, Y] = meshgrid(1:selectedK,1:selectedK);
    % Display the background
    imagesc(I);
    hold on;
    % Set the background color to white
    set(gcf,'Colormap',[1, 1, 1]);
    % Insert the labels
    for c_out = 1:selectedK
        for c_in = 1:selectedK
            if TM_pval2sided(s_ind,c_out,c_in) < 0.05 && TM_pval2sided(s_ind,c_out,c_in) > (0.05/selectedK)
                text(X(c_out,c_in),Y(c_out,c_in)+.22,'*','FontSize',16,'HorizontalAlignment','center','FontWeight','bold','Color','r');
            end
            if TM_pval2sided(s_ind,c_out,c_in) < (0.05/selectedK) && TM_pval2sided(s_ind,c_out,c_in) > (0.05/(selectedK*selectedK))
                text(X(c_out,c_in),Y(c_out,c_in)+.22,'*','FontSize',16,'HorizontalAlignment','center','FontWeight','bold','Color','g');
            end
            if TM_pval2sided(s_ind,c_out,c_in) <= (0.05/(selectedK*selectedK))
                text(X(c_out,c_in),Y(c_out,c_in)+.22,'*','FontSize', 16,'HorizontalAlignment','center','FontWeight','bold','Color','b');
            end
        end
    end
    % Calculte the grid lines
    grid = 0.5:1:selectedK+0.5;
    grid1 = [grid; grid];
    grid2 = repmat([0.5;selectedK+0.5],1,length(grid));
    % Plot the grid lines
    plot(grid1,grid2,'k')
    plot(grid2,grid1,'k')
    ticklabels = cell(1,selectedK);
    for p = 1:selectedK
        ticklabels{p} = p;
    end
    xticks(1:1:selectedK); xticklabels(ticklabels);
    yticks(1:1:selectedK); yticklabels(ticklabels);
    set(gca,'FontSize',9)
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel('From FC State','FontSize',11)
    xlabel('To FC State','FontSize',11)
    axis square
    box off
    
    sgtitle([{'Two-sided {\itp}-value'},{'State-to-state transition probability'}],'Fontsize',14,'FontWeight','bold') 
end
saveas(Fig, fullfile(data_dir, ['K' num2str(selectedK) '_TransitionsDifferences.png']),'png');
saveas(Fig, fullfile(data_dir, ['K' num2str(selectedK) '_TransitionsDifferences.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_TransitionsDifferences']);

close all;