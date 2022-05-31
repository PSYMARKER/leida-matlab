function Plot_K_nodes_in_cortex(data_dir,save_dir,selectedK,parcellation)
%
% Plot the nodes with positive values in the vectors of the centroids red
% and the remaning nodes blue. Scale the size and transparency of the nodes
% according to their contribution.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot the nodes according to their sign and contribution for
%               each centroid
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
V = V./max(abs(V));
clear Kmeans_results

% Get centers of gravity of each brain area from the parcellation provided
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
sz = size(Volume);
MNI_coord = zeros(max(Volume(:)),3);
for n = 1:max(Volume(:))
    ind_N = find(Volume == n);
    [X,Y,Z] = ind2sub(sz,ind_N);
    MNI_coord(n,:) = mean([Y,X,Z]);
end
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

disp(' ');
disp(['Plotting nodes in cortex for the ' num2str(selectedK) ' PL states:'])
Fig = figure('Position', get(0, 'Screensize'));

for c = 1:selectedK
    % Select row c of the set of centroids
    Vc = V(c,:);
    
    subplot_tight(1,selectedK,c,0.01)
    hold on

    % Cortex patch
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
    reducepatch(cortexpatch,0.5);
    isonormals(scortex,cortexpatch);

    a = 2;
    [x,y,z] = sphere;
    x = a*x;
    y = a*y;
    z = a*z;

    for n = 1:length(Vc)
        if Vc(n) > 0 % color nodes with positive contribution in a gradient of red colors
            surf(Vc(n)*x+MNI_coord(n,1), Vc(n)*y+MNI_coord(n,2), Vc(n)*z+MNI_coord(n,3),'FaceColor',[Vc(n) 0 0],...,
                'EdgeColor','none','FaceAlpha',abs(Vc(n)));
        elseif Vc(n) < 0 % color nodes with negative contribution in a gradient of blue colors
            surf(Vc(n)*x+MNI_coord(n,1), Vc(n)*y+MNI_coord(n,2), Vc(n)*z+MNI_coord(n,3),'FaceColor',[0 0 -Vc(n)],...,
                'EdgeColor','none','FaceAlpha',abs(Vc(n)));
        end
    end
    
    if c == 1
        title('Global Mode', 'Fontsize', 12);
    else
        title(['PL state ' num2str(c)], 'Fontsize', 12)
    end

    axis off;
    axis equal;

    material dull;
    lighting gouraud;
    daspect([1 1 1])
    view(-90,90)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    
end

saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_NodesCortex.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_NodesCortex.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_NodesCortex']);

close all;