function Plot_K_matrix(data_dir,save_dir,selectedK)
%
% Plot the centroids in matrix format
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
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
clear Kmeans_results

% Get the number of areas based on centroids size
n_areas = size(V,1);

% Reorder the position of each parcel
Order = [1:2:n_areas n_areas:-2:2];

disp(' ');
disp(['Plotting the ' num2str(selectedK) ' PL states in matrix format:'])
Fig = figure('Position', get(0, 'Screensize'));  
for c = 1:selectedK
    subplot_tight(1,selectedK,c,0.02)
    colormap(jet)
    imagesc(V(Order,c)*V(Order,c)',[-1 1])
    % set(gca, 'XTick',10:20:n_areas, 'YTick',10:20:n_areas)
    axis square
    title(['V_{C_{' num2str(c) '}}' '.V_{C_{' num2str(c) '}}' '^T'], 'Fontsize', 12)
end
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_Matrix.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_Matrix.fig']),'fig');
disp(['- Plot of matrices successfully saved as K' num2str(selectedK) '_Matrix']);

close all;