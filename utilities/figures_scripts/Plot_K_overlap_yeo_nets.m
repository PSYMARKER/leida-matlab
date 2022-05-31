function Plot_K_overlap_yeo_nets(data_dir,save_dir,selectedK,parcellation)
%
% Compute the overlap between each PL state (1 to selectedK) and the RSNs
% defined in Yeo et al., (2011). Plot the overlap between each PL state and
% each of the Yeo RSNs.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of the overlap of each PL state with each of the RSNs
%               defined in Yeo et al., (2011)
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Matrix of dimension selectedK*90 where each row represents one FC state
V = Kmeans_results{rangeK == selectedK}.C;

% Number of areas considered for analysis
n_areas = size(V,2);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
% Name of the Yeo et al., (2011) networks
YeoNames = {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

disp(' ');
disp(['Plotting overlap with Yeo RSNs for K = ' num2str(selectedK) ' clusters:'])
Fig = figure('Position', get(0, 'Screensize'));

b = bar(1:selectedK,squeeze(cc_V_yeo7(rangeK == selectedK,1:selectedK,1:7)));
for i = 1:7
    set(b(i),'FaceColor',YeoColor(i,:));
end
xlabel(['K = ' num2str(selectedK) ' clusters'], 'FontSize',14);
ylabel({'Pearson''s correlation with 7 RSNs';'(Yeo et al. 2011)'}, 'FontSize', 14);
legend('Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode', ...
       'Location','northwest','FontSize',14);
legend boxoff;
box off;
for centroidK = 1:selectedK
    for NetYeo = 1:7
        if p_V_yeo7(rangeK == selectedK,centroidK,NetYeo) < 0.05/selectedK && ~isnan(p_V_yeo7(rangeK == selectedK,centroidK,NetYeo))
            if isprop(b(1), 'XEndPoints')
                x_K = b(NetYeo).XEndPoints;
                x_in_group = x_K(centroidK);
                y_K = b(NetYeo).YEndPoints;
                y_in_group = y_K(centroidK);
                if cc_V_yeo7(rangeK == selectedK,centroidK,NetYeo) > 0
                    text(x_in_group,y_in_group*1.01,'*','fontsize',25,'Color',[0 0 0],'HorizontalAlignment','center');
                else
                    text(x_in_group,0.01,'*','fontsize',25,'Color',[0 0 0],'HorizontalAlignment','center');
                end
            else
                % If in this Matlab version it was not possible to add asterisks on the bars
                % These can be added manually at the following locations:
                disp(['- Centroid ' num2str(centroidK) ' significantly overlaps with ' YeoNames{NetYeo} ' Yeo RSN:'])
                disp(['     cc = ' num2str(cc_V_yeo7(rangeK == selectedK,centroidK,NetYeo)) ' p = ' num2str(p_V_yeo7(rangeK == selectedK,centroidK,NetYeo))])
            end
        end
    end
end           
     
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_OverlapYeoNets.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_OverlapYeoNets.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_OverlapYeoNets']);

close all;