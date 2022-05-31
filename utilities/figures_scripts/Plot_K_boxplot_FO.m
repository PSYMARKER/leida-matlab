function Plot_K_boxplot_FO(data_dir,save_dir,selectedK,parcellation)
%
% Plot the boxplot of fractional occupancy (FO) values by condition for
% each of the centroids given by the specified K.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation atlas used to segement the brain
%
% OUTPUT:
% .fig/.png     Boxplot of the fractional occupancy values for each
%               condition across centroids
%
% Authors: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end
if isfile([data_dir file_P])
    load([data_dir file_P], 'cond', 'P', 'P_pval2sided', 'Index_Conditions');
end

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Matrix of dimension selectedK*n_areas where each row represents one PL state
V = Kmeans_results{rangeK == selectedK}.C;

% Number of areas considered for analysis
n_areas = size(V,2);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);

disp(' ');
disp(['Plotting the boxplot of fractional occupancy values for each of the ' num2str(selectedK) ' PL states by condition:'])
Fig = figure('Position', get(0, 'Screensize'));
for c = 1:selectedK
    
    if selectedK <= 10
        subplot(4,selectedK,[c+selectedK c+2*selectedK])
    else
        subplot(2,10,c)
    end
    
    [~, net] = max(cc_V_yeo7(rangeK == selectedK,c,:));
    P_cond = cell(1,n_Cond);
    g = cell(1,n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j,rangeK == selectedK,c);
        g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == selectedK,c)),1);
    end
    P_data = vertcat(P_cond{:});
    g_data = vertcat(g{:});
    
    bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
    set(bx,{'linew'},{1})
    b = get(get(gca,'children'),'children');   % Get the handles of all the objects
    f = get(b,'tag');   % List the names of all the objects
    med_col = b(n_Cond+1:n_Cond*2);
    set(med_col, 'Color', 'k');
    box_col = b(n_Cond*2+1:n_Cond*3);
    if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
        set(box_col, 'Color', YeoColor(net,:));
    else
        set(box_col, 'Color', 'k');
    end   
    
    hold on
    set(gca,'XTickLabel',cond)
    xtickangle(30)
    if c == 1 || c == 11
        ylabel('Fractional Occupancy (%)', 'Fontsize', 10);
    end
    title(['PL State ' num2str(c)]);
    set(gca,'color','none')
    
    X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
    cond_pair = 1;
    for cond1 = 1:n_Cond-1
        for cond2 = cond1+1:n_Cond
            X_locations(cond_pair,1) = cond1;
            X_locations(cond_pair,2) = cond2;
            cond_pair = cond_pair + 1;
        end
    end
    % X_locations = sortrows(X_locations,2);
    X_locations(:,1) = X_locations(:,1) + 0.1;
    X_locations(:,2) = X_locations(:,2) - 0.1;
    Max_Y =  max(P_data);
    Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(P_pval2sided(:,rangeK == selectedK,c) <= 0.05/selectedK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == selectedK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(P_pval2sided(:,rangeK == selectedK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(P_pval2sided(asterisks,rangeK == selectedK,c),'%10.1e')],'Fontsize',7)
    % end
    
    ylim([0 max(Y_LOCATIONS)+0.15])
    
    hold off
    box off
end
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_BoxplotFracOccup.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_BoxplotFracOccup.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_BoxplotFracOccup']);

close all;
