function Plot_K_V1_VoxelSpace(leida_dir,selectedK,parcellation)
%
% Plot the mean centroids computed in MNI 10mm voxel space.
%
% INPUT:
% leida_dir     directory where LEiDA results are stored
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of centroids rendered in 3D glass brain in MNI 10mm
%               voxel space colored according to overlap with Yeo RSNs
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% INPUT EXAMPLES:
% leida_dir = 'D:/LEiDA_Toolbox/LEiDA_Results_ABIDE_dparsf_AAL120/';
% selectedK = 15;
% parcellation = 'AAL120';

% Check whether a folder with the results for selectedK exists
if ~exist([leida_dir 'K' num2str(selectedK) '/'], 'dir')
    mkdir([leida_dir 'K' num2str(selectedK) '/']);
end
K_dir = [leida_dir 'K' num2str(selectedK) '/'];
% Open the directory where the results are saved
cd(K_dir);

% File with the results from voxel space analysis (output from EigenVectors_VoxelSpace.m)
file_V1_VoxelSpace = 'V1_VoxelSpace.mat';

% Load required data:
if isfile(file_V1_VoxelSpace)
    load(file_V1_VoxelSpace, 'mean_V1', 'ind_voxels', 'MNI10mm_Mask');
end

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([leida_dir file_cluster])
    load([leida_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Matrix (selectedK*n_areas) where each row represents one PL state
V = Kmeans_results{rangeK == selectedK}.C;
% Get the number of areas based on centroids size
n_areas = size(V,2);

% Load 2mm MNI voxel space (to define a transparent brain)
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Color code from paper Yeo et al., 2011
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Patch_thr=0.5;

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

disp(' ');
disp(['Plotting the ' num2str(selectedK) ' PL states in voxel space:'])
Fig = figure('Position', get(0, 'Screensize'));

for c = 1:selectedK
    
    % Get the net with which PL state overlaps the most
    [~, net] = max(cc_V_yeo7(rangeK == selectedK,c,:));
    
    % Get the vector representing this centroid
    V_state = zeros(size(MNI10mm_Mask));
    V_state(ind_voxels) = mean_V1(c,:);
    V_state = imresize3(V_state, [91 109 91]);
    V_state(~scortex) = 0;
    
    subplot_tight(3,K,c,0.02)
    
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);
    
    hold on
    if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', YeoColor(net,:),...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    else
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', 'k',...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    reducepatch(Centroid_patch, 0.1);
    isonormals(V_state, Centroid_patch);
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight
    xlim([5 105])
    ylim([5 85])
    zlim([0 80])
    axis off
    
    subplot_tight(3,K,c+K,0.02)
    
    cortexpatch = patch(isosurface(scortex,0.5), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);
    
    hold on
    if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', YeoColor(net,:),...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    else
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', 'k',...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    isonormals(V_state,Centroid_patch);
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)    % Top view    Side view:   view(0,0)
    xlim([5 105])
    ylim([5 85])
    zlim([0 80])
    camlight
    axis off
    
    subplot_tight(3,K,c+2*K,0.02)
    
    cortexpatch=patch(isosurface(scortex,0.5), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.1);
    isonormals(scortex,cortexpatch);
    
    hold on
    if p_V_yeo7(rangeK == selectedK,c,net) < 0.05/selectedK
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', YeoColor(net,:),...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    else
        Centroid_patch = patch(isosurface(V_state, Patch_thr*max(mean_V1(c,:))), 'FaceColor', 'k',...,
                               'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    isonormals(V_state,Centroid_patch);
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(90,0)    % Top view    Side view:   view(0,0)
    xlim([5 105])
    ylim([5 85])
    zlim([0 80])
    camlight
    axis off
    
end

saveas(Fig, fullfile(K_dir, ['K' num2str(selectedK) '_V1_VoxelSpace.png']),'png');
saveas(Fig, fullfile(K_dir, ['K' num2str(selectedK) '_V1_VoxelSpace.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_V1_VoxelSpace']);

close all;