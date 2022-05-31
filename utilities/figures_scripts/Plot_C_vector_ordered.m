function Plot_C_vector_ordered(data_dir,save_dir,selectedK,centroid,parcellation)
%
% Plot the selected centroid in vector format with the labels of each area
% and ordered by contribution
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% centroid      centroid defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot centtroid as a barplot and label each area
%               according to the applied parcellation and order areas by
%               contribution
%
% Authors: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

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

% Get labels from parcellation used
switch parcellation
    case 'AAL116'
        load('ParcelsMNI2mm','label116');
        labels = label116([1:n_areas],:);
        clear label116
    case 'AAL120'
        load('ParcelsMNI2mm','label120');
        labels = label120([1:n_areas],:);
        clear label120    
end

[Vo, asc_ind] = sort(V(:,centroid));
numMarkers = size(Vo,1);
if centroid == 1
    markerColors = jet(numMarkers*2);
else
    markerColors = jet(numMarkers);
end
miny = min(Vo);
maxy = max(Vo);
colorMapRows = round(rescale((Vo - miny) / (maxy - miny), 1, numMarkers));

disp(' ');
disp(['Plotting the vector format with areas ordered by contribution for the PL state ' num2str(centroid) ':'])
Fig = figure('Position', get(0, 'Screensize'));
subplot(1,5,3)
for k = 1:length(Vo)
    thisMarkerColor = markerColors(colorMapRows(k), :);
    barh(k,Vo(k),'FaceColor',thisMarkerColor,'EdgeColor','none','Barwidth',.5)
    hold on;
    ylim([0 n_areas+1])
    xlim([-1 1])
    % set(gca,'XTick',[-1 0 1])
    set(gca,'YTick',1:n_areas)
    set(gca,'YTickLabel',labels(asc_ind,:),'Fontsize',6)
    ax = gca;
    ax.XAxis.FontSize = 10;
    grid on;
    title(['PL state ' num2str(centroid)],'Fontsize',12)
end

saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_VectorOrdered.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_VectorOrdered.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) 'C' num2str(centroid) '_VectorOrdered']);

close all;