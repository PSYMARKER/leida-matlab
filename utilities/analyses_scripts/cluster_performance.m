function cluster_performance(data_dir)
%
% Compute the Dunn's index, Calinski-Harabasz (CH) index and average 
% Silhouette coefficient for each value of K to assess clustering
% performance.
%
% INPUT:
% data_dir      directory where the results from K-means are saved
%
% OUTPUT:
% dunn_score    Dunn's index computed for each K
% avg_sil       average Silhouette coefficient computed for each K
% CH            CH index computed for each K
% .fig/.png     plot of the Dunn's index, average Silhouette
%               coefficient and CH index for each number of clusters
%
% Author: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%         Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% Input example:
% data_dir = 'D:/LEiDA_Toolbox/LEiDA_results/';

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_V1])
    load([data_dir file_V1], 'V1_all');
end
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Dunn's index
disp(' ')
disp('Computing Dunn''s index:')
distM = squareform(pdist(V1_all,'cosine'));
dunn_score = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    dunn_score(i) = dunns(rangeK(i), distM, Kmeans_results{i}.IDX);
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxdunn] = max(dunn_score);
disp(['- Best clustering solution according to Dunn''s index: ' num2str(rangeK(ind_maxdunn)) ' clusters']);

% Average Silhouette Coefficient
disp(' ')
disp('Computing average Silhouette coefficient:')
avg_sil = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_sil = evalclusters(V1_all,Kmeans_results{i}.IDX','Silhouette','Distance','cosine');
    avg_sil(i) = eva_sil.CriterionValues;
    clear eva_sil;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxsil] = max(avg_sil);
disp(['- Best clustering solution according to average Silhouette coefficient: ' num2str(rangeK(ind_maxsil)) ' clusters']);

% CH index
disp(' ')
disp('Computing CH index:')
CH = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_CH = evalclusters(V1_all,Kmeans_results{i}.IDX','CalinskiHarabasz');
    CH(i) = eva_CH.CriterionValues;
    clear eva_CH;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxCH] = max(CH);
disp(['- Best clustering solution according to CH index: ' num2str(rangeK(ind_maxCH)) ' clusters']);

% Saving results from clustering performance analysis
save_file = 'ClusterPerformance.mat';
save([data_dir '/' save_file],'CH','dunn_score','avg_sil');
disp(' ')
disp(['Clustering performance results saved successfully as ' save_file]);
disp(' ')

disp('Plotting clustering performance results:')
Fig = figure('Position', get(0, 'Screensize'));
x = 2:1:rangeK(end);
tiledlayout(3,1);
ax1 = nexttile;
plot(ax1,x,dunn_score,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax1,'Dunn''s index','Fontsize',12);
box off;

ax2 = nexttile;
plot(ax2,x,avg_sil,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax2,'Average Silhouette coefficient','Fontsize',12);
box off;

ax3 = nexttile;
plot(ax3,x,avg_sil,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax3,'CH index','Fontsize',12);
box off;
xlabel('Number of clusters','Fontsize',12);

saveas(Fig, fullfile(data_dir, 'ClusterPerformance.png'),'png');
saveas(Fig, fullfile(data_dir, 'ClusterPerformance.fig'),'fig');
disp('- Plot successfully saved as ClusterPerformance');

close all;
