function cluster_stability(data_dir,save_dir,selectedK)
%
% Compute the clustering stability for the selected K using
% the procedure explained in Farinha et al., (2022).
%
% INPUT:
% data_dir      directory where the results from K-means are saved
% save_dir      directory to save the results, it should correspond
%               to the folder used to store results for selectedK
% selectedK     number of clusters
%
% OUTPUT:
% ConfMat_fold  confusion matrix for each fold
% acc           accuracy computed from the confusion matrix
% VI_fold       variation of information
% rand_fold     adjusted rand index
% figure        tables with cluster stability results
%
% Author: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
% Input examples:
% data_dir = 'D:/LEiDA_Toolbox/LEiDA_Results_ABIDE_dparsf_AAL120/';
% save_dir = 'D:/LEiDA_Toolbox/LEiDA_Results_ABIDE_dparsf_AAL120/K15/';
% selectedK = 15;

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
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING STABILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Create K-folds data for K-fold cross-validation (10 folds)
rng('default') % For reproducibility
cv = cvpartition(size(V1_all,1),'KFold',10);

% Agora aplicamos o procedimento para cada fold
numFolds = cv.NumTestSets; % number of folds
ConfMat_fold = cell(1,numFolds); % confusion matrix of each fold
acc = zeros(1,numFolds); % accuracy for each fold
VI_fold = zeros(1,numFolds); % variation of information for each fold
rand_fold = zeros(1,numFolds); % adjusted rand index for each fold
Kmeans_results = cell(1,numFolds); % K-means results for each fold

disp(' ')
disp(['Assessing stability of ' num2str(selectedK) ' clusters using ' num2str(numFolds) '-fold cross-validation:'])
for i = 1:numFolds
    disp(['- Fold ' num2str(i)]);
    idx_train = cv.training(i);
    trainData = V1_all(idx_train,:);
    
    % Apply K-means to the training data of fold i -> Partition P1
    disp('      - Applying K-means to training dataset');
    [train_IDX,train_C,~,~] = kmeans(trainData,selectedK,'Distance','cosine','Replicates',500,'MaxIter',1500,...,
                                   'OnlinePhase','off','Display','off','Options',statset('UseParallel',0));
    
    % Relabel cluster numbers according to probability of occurrence
    [~, ind_sort] = sort(hist(train_IDX,1:selectedK),'descend');
    train_C = train_C(ind_sort,:);
    Kmeans_results{i}.Ctrain = train_C; 
    
    % Get labels of test data using the nearest centroid classifier -> P2
    idx_test = cv.test(i);
    testData = V1_all(idx_test,:);
    disp('      - Labelling test dataset using nearest centroid classifier');
    for j = 1:size(testData)
        [~,label_test] = pdist2(train_C,testData,'cosine','Smallest',1);
    end
    Kmeans_results{i}.label = label_test;
    
    % Apply K-means to the test data of fold i -> Partition P3
    disp('      - Applying K-means to test dataset');
    [test_IDX,~,~,~] = kmeans(testData,selectedK,'Distance','cosine','Replicates',500,'MaxIter',1500,...,
                            'OnlinePhase','off','Display','off','Options',statset('UseParallel',0));
    
    % Relabel cluster numbers according to probability of occurrence
    [~, ind_sort] = sort(hist(test_IDX,1:selectedK),'descend');
    [~,idx_sort] = sort(ind_sort,'ascend');
    test_IDX = idx_sort(test_IDX);
    Kmeans_results{i}.IDXtest = test_IDX;
    
    % Compute confusion matrix for each fold
    % P2 -> partition from classifier: "known" groups
    % P3 -> partition from clustering of test sample: "predicted" groups
    disp('      - Computing confusion matrix');
    ConfMat = confusionmat(label_test,test_IDX);
    ConfMat_fold{i} = ConfMat;
    acc(:,i) = trace(ConfMat) / sum(ConfMat(:));
    disp(['      - Accuracy: ' num2str(acc(i))]);
    
    % Variation of information
    n = length(label_test);
    sigma = 0.0;
    for m = 1:selectedK
        p = length(find(label_test == m)) / n;
        for a = 1:selectedK
            q = length(find(test_IDX == a)) / n;
            r = length(intersect(find(label_test == m),find(test_IDX == a))) / n;
            if r > 0.0
                sigma = sigma + r * (log(r / p) + log(r / q));
            end
        end
    end
    VI_fold(:,i) = abs(sigma);
    disp(['      - Variation of information: ' num2str(VI_fold(i))]);
    
    rand_fold(:,i) = rand_index(label_test, test_IDX, 'adjusted');
    disp(['      - Adjusted Rand Index: ' num2str(rand_fold(i))]);
    disp(' ')
end

% Saving results from clustering performance analysis
save_file = ['K' num2str(selectedK) '_ClusterStability.mat'];
save([save_dir '/' save_file],'ConfMat_fold','acc','VI_fold','rand_fold','Kmeans_results','cv');
disp(['Clustering stability results for ' num2str(selectedK) ' clusters saved successfully as ' save_file]);

% Plot accuracy, variation of information and adjusted rand index
disp(' ')
disp('Plotting the accuracy, variation of information and adjusted rand index clustering stability results:')
Fig = figure('Position', get(0, 'Screensize'));
clust_stab = [acc; VI_fold; rand_fold];
imagesc(clust_stab, 'AlphaData', 0.8);
hold on
set(gcf,'Colormap',[1, 1, 1]);
for c = 1:size(clust_stab,1)
    for j = 1:numFolds
        caption = sprintf('%.3f', clust_stab(c,j));
        text(j,c,caption,'FontSize',15,'Color',[0, 0, 0],'FontWeight','bold','HorizontalAlignment','Center');
    end
end
grid_h = 0.5:1:3+0.5;
gridh1 = [grid_h; grid_h];
gridh2 = repmat([0.5;numFolds+0.5],1,length(grid_h));
% Plot the grid lines
plot(gridh2,gridh1,'k')
grid_v = 0.5:1:numFolds+0.5;
gridv1 = [grid_v; grid_v];
gridv2 = repmat([0.5;3+0.5],1,length(grid_v));
plot(gridv1,gridv2,'k')
ticklabels = cell(1,numFolds);
for p = 1:numFolds
    ticklabels{p} = p;
end
yticks(1:1:3); yticklabels({'Accuracy', 'Variation of information', 'Adjusted rand index'});
xticks(1:1:numFolds); xticklabels(ticklabels);
set(gca,'FontSize',14,'FontWeight','bold')
axis square
box off
  
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterStability.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterStability.fig']),'fig');
disp(['- Plot of the clustering stability results successfully saved as K' num2str(selectedK) '_ClusterStability']);

close all;