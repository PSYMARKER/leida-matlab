function Plot_K_repertoire(data_dir,save_dir,selectedK,parcellation)
%
% Plot the centroids in 3D glass brain, in vector format and the boxplots
% of the fractional occupancy and dwell time value.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of summary information for the repertoire of K
%               centroids
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';
% File with results for the dwell time (output from LEiDA_stats_DwellTime.m)
file_LT = 'LEiDA_Stats_DwellTime.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end
if isfile([data_dir file_P])
    load([data_dir file_P], 'cond', 'P', 'P_pval2sided', 'Index_Conditions');
end
if isfile([data_dir file_LT])
    load([data_dir file_LT], 'LT', 'LT_pval2sided');
end

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Matrix of dimension selectedK*90 where each row represents one FC state
V = Kmeans_results{rangeK == selectedK}.C;
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));

% Number of areas considered for analysis
n_areas = size(V,1);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Yeo_names = {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Reorder the position of each parcel
Order = [1:2:n_areas n_areas:-2:2];

disp(' ');
disp(['Plotting summary information for the ' num2str(selectedK) ' PL states:'])
Fig = figure('Position', get(0, 'Screensize')); 
colormap(jet)
for c = 1:selectedK
    
    [~, net] = max(cc_V_yeo7(rangeK == selectedK,c,:));
    
    % First Plot view from top
    subplot_tight(9,selectedK,c+selectedK,0.02)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to contribution
    n_pos = find(V(:,c) > 0);
    if numel(n_pos) > 0
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else % Plot of the global mode with "gray/transparent" color
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
        title({['FC state ' num2str(c), newline, Yeo_names{net}]} ) %,Yeo_names{c}})
    elseif numel(n_pos) == 0
        title({['FC state ' num2str(c), newline, 'Global Mode']} ) %,Yeo_names{c}})
    else
        title({['FC state ' num2str(c), newline, '']} ) %,Yeo_names{c}})
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    subplot_tight(9,selectedK,c+2*selectedK,0.025)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    n_pos = find(V(:,c) > 0);
    if n_pos
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol = Volume > 0 & Volume <= 90; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    % PLOT THE VECTOR OF EACH FC STATE
    subplot_tight(9,selectedK,[c+3*selectedK c+4*selectedK],0.045)
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'XTick',[-1,0,1],'Fontsize',8)
    set(gca,'YTick',[10:20:n_areas],'Fontsize',8)
    ax = gca;
    ax.XAxis.FontSize = 8;
    grid on
    
    % BOXPLOT FRACTIONAL OCCUPANCY
    subplot_tight(9,selectedK,[c+5*selectedK c+6*selectedK],0.035)
    
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
    set(gca,'XTickLabel',cond, 'Fontsize',6)
    xtickangle(30)
    if c == 1
        ylabel('Fractional Occupancy (%)');
    end
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
    
    % BOXPLOT DWELL TIME
    subplot_tight(9,selectedK,[c+7*selectedK c+8*selectedK],0.035)
    
    LT_cond = cell(1,n_Cond);
    r = cell(1,n_Cond);
    for j = 1:n_Cond
        LT_cond{j} = LT(Index_Conditions == j,rangeK == selectedK,c);
        r{j} = repmat(cond(j), length(LT(Index_Conditions == j,rangeK == selectedK,c)),1);
    end
    LT_data = vertcat(LT_cond{:});
    r_data = vertcat(r{:});
    
    bx = boxplot(LT_data, r_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
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
    set(gca,'XTickLabel',cond, 'Fontsize',6)
    xtickangle(30)
    if c == 1
        ylabel('Dwell Time (s)');
    end
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
    Max_Y =  max(LT_data);
    Y_LOCATIONS = Max_Y + (abs(X_locations(:,1) - X_locations(:,2)));
    
    % Green asterisks
    asterisks = find(LT_pval2sided(:,rangeK == selectedK,c) <= 0.05/selectedK);
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.035, '*','Color','g','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
    % if asterisks
        % text(mean(X_locations(asterisks,:),2) + 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
            % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == selectedK,c),'%10.1e')],'Fontsize',7)
    % end
    
    % Blue asterisks
    asterisks = find(LT_pval2sided(:,rangeK == selectedK,c) <= 0.05/sum(rangeK));
    plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
    plot(mean(X_locations(asterisks,:),2), Y_LOCATIONS(asterisks)*1.035, '*b','Markersize',4)
    
    % Uncomment to add p-values to plot (add -0.05 in the line above)
     % if asterisks
         % text(mean(X_locations(asterisks,:),2) - 0.05 , Y_LOCATIONS(asterisks)*1.05,...,
             % [repmat('p=',numel(asterisks),1) num2str(LT_pval2sided(asterisks,rangeK == selectedK,c),'%10.1e')],'Fontsize',7)
     % end
    
    ylim([0 max(Y_LOCATIONS)*1.15])
    
    hold off
    box off
end
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_RepertoireSummary.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_RepertoireSummary.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_RepertoireSummary']);

close all;   
