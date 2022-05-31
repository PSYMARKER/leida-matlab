function [Kmeans_results,rangeK] = LEiDA_cluster(data_dir)
%
% Cluster all leading eigenvectors into a set of K clusters using
% the K-means algorithm. 
% This function returns an optimal solution for each values of K.
%
% INPUT:
% data_dir         directory where the LEiDA_EigenVectors file is saved
%
% OUTPUT:
% Kmeans_results   clustering solutions for the range of values K
% rangeK           range of values of K considered for clustering
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% Set maximum/minimum number of clusters
mink = 2;
maxk = 20;

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';

% Load the leading eigenvectors
load([data_dir file_V1], 'V1_all');

% Range of clustering solutions to be computed
rangeK = mink : maxk;

% Data type to store results for each K using the K-means algorithm
Kmeans_results = cell(size(rangeK));

% Number of new initial cluster centroid positions to run
% (it is convenient to run more replicates for larger samples)
replicates = ceil(50 + size(V1_all, 1)/300);

disp(' ');
disp('Clustering eigenvectors into:')

% For each K do:
for K = 1:length(rangeK)
    disp(['- ' num2str(rangeK(K)) ' FC states'])
    [IDX, C, SUMD, D] = kmeans(V1_all,rangeK(K),'Distance','Cosine','Replicates',replicates,'MaxIter',1000,'OnlinePhase','off',...,
                              'Display','off','Options',statset('UseParallel',1));
                         
    % ind_sort sorts the clusters in descending order of occupancy
    [~, ind_sort] = sort(hist(IDX,1:rangeK(K)),'descend');
    [~, idx_sort] = sort(ind_sort,'ascend');
    Kmeans_results{K}.IDX = idx_sort(IDX);     % Cluster time course - numeric collumn vector
    Kmeans_results{K}.C = C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results{K}.SUMD = SUMD(ind_sort);   % Within-cluster sums of point-to-centroid distances
    Kmeans_results{K}.D = D(:,ind_sort);       % Distance from each point to every centroid
end

% Name of the file to save output
save_file = 'LEiDA_Clusters.mat';

save([data_dir '/' save_file], 'Kmeans_results', 'rangeK')
disp(['K-means clustering completed and results saved as ' save_file])
disp(' ');