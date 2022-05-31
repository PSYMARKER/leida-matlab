function Plot_K_3Dbrain(data_dir,save_dir,selectedK,parcellation)
%
% Plot the centroids in a 3D glass brain indicating their overlap with the
% networks defined by Yeo et al., (2011).
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of centroids rendered in 3D glass brain stating
%               the overlap of each with the Yeo networks
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Matrix (selectedK*n_areas) where each row represents one PL state
V = Kmeans_results{rangeK == selectedK}.C;
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));

% Get the number of areas based on centroids size
n_areas = size(V,1);

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Yeo_names = {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

disp(' ');
disp(['Plotting the ' num2str(selectedK) ' PL states in 3D glass brain:'])
Fig = figure('Position', get(0, 'Screensize'));

for c = 1:selectedK
    
    [~, net] = max(cc_V_yeo7(rangeK == selectedK,c,:));
    
    % First Plot view from top
    if selectedK <= 10
        subplot_tight(2,selectedK,c,0.01)
    else
        if c <= 10
            subplot_tight(4,10,c,0.03)
        else
            subplot_tight(4,10,10+c,0.03)
        end
    end
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to
    % contribution
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
        title({['PL state ' num2str(c), newline, Yeo_names{net}]}, 'Fontsize', 12) %,Yeo_names{c}})
    elseif numel(n_pos) == 0
        title({['PL state ' num2str(c), newline, 'Global Mode']}, 'Fontsize', 12) %,Yeo_names{c}})
    else
        title({['PL state ' num2str(c), newline, '']}, 'Fontsize', 12) %,Yeo_names{c}})
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    if selectedK <= 10
        subplot_tight(2,selectedK,c+selectedK,0.01)
    else
        if c <= 10
            subplot_tight(4,10,10+c,0.02)
        else
            subplot_tight(4,10,20+c,0.02)
        end
    end
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
    
end

saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_3Dbrain.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_3Dbrain.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_3Dbrain']);

close all;