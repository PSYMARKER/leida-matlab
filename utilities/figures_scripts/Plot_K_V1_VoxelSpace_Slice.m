function Plot_K_V1_VoxelSpace_Slice(leida_dir,selectedK)
%
% Plot the mean centroids computed in MNI 10mm voxel space at each slice.
%
% INPUT:
% leida_dir     directory where LEiDA results are stored
% selectedK     K defined by the user
%
% OUTPUT:
% .fig/.png     Plot of centroids in different slices in MNI 10mm voxel
%               space
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% INPUT EXAMPLES:
% leida_dir = 'D:/LEiDA_Toolbox/LEiDA_Results_ABIDE_dparsf_AAL120/';
% selectedK = 15;

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

% On slices
plane = 'horizontal';

disp(' ');
disp(['Plotting the ' num2str(selectedK) ' PL states in slice (' plane '):'])
figure
colormap(jet)

Slices=9:-1:2;

State_Volumes = zeros(selectedK,size(MNI10mm_Mask,1)*3,size(MNI10mm_Mask,2)*3,size(MNI10mm_Mask,3)*3);

c = 1;
for Cluster = 2:selectedK
    
    V_state = zeros(size(MNI10mm_Mask));   
    V_state(ind_voxels) = mean_V1(Cluster,:);
    V_state = imresize3(V_state,3);
    
    State_Volumes(Cluster,:,:,:) = V_state;
    lim = max(abs(mean_V1(Cluster,:)));
    
    %V_state(Brain_Mask_low==0)=-1;
    
    for z = 1:length(Slices) - 1
        
        subplot(length(Slices)-1,selectedK-1,c+(z-1)*(selectedK-1))
        
        if strcmp(plane,'horizontal')
            slice_to_plot = round(size(V_state,3)/(length(Slices)+2))*Slices(z+1);
            imagesc(squeeze(V_state(:,:,slice_to_plot))',[-lim lim])
        elseif strcmp(plane,'vertical')
            slice_to_plot = round(size(V_state,2)/(length(Slices)+2))*Slices(z+1);
            imagesc(squeeze(V_state(:,slice_to_plot,:))',[-lim lim])
        elseif strcmp(plane,'side')
            slice_to_plot = round(size(V_state,1)/(length(Slices)+2))*Slices(z+1);
            imagesc(squeeze(V_state(slice_to_plot,:,:))',[-lim lim])
        end
        
        axis equal
        axis tight
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        axis xy
        if z == 1
            title(['V_{C' num2str(Cluster) '}']);
        end
        
        if Cluster == 2
            ylabel(['Z = ' num2str(Slices(z))])
        end
        view(0,90)
    end
    c = c+1;
end

saveas(Fig, fullfile(K_dir, ['K' num2str(selectedK) '_V1_VoxelSlice.png']),'png');
saveas(Fig, fullfile(K_dir, ['K' num2str(selectedK) '_V1_VoxelSlice.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_V1_VoxelSlice']);

close all;

