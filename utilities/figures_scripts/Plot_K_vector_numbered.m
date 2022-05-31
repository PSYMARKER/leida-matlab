function Plot_K_vector_numbered(data_dir,save_dir,selectedK)
%
% Plot the centroids in vector format with the number of each area.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
%
% OUTPUT:
% .fig/.png     Plot each centroid as a barplot and number of each area
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
disp(['Plotting the ' num2str(selectedK) ' PL states as barplot with numbered areas:'])
Fig = figure('Position', get(0, 'Screensize')); 
for c = 1:selectedK
    if selectedK <= 10
        subplot(1,selectedK,c)
    else
        subplot_tight(2,10,c,0.065)
    end
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo < 0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo >= 0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'YTick',[10:10:n_areas], 'Fontsize',12)
    ax = gca;
    ax.XAxis.FontSize = 12;
    grid on
    title(['PL state ' num2str(c)], 'Fontsize',12)
end
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_VectorNumber.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_VectorNumber.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_VectorNumber']);

close all;