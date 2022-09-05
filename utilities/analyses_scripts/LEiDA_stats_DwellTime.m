function LEiDA_stats_DwellTime(data_dir,cond,pair,tr,n_permutations,n_bootstraps)
%
% For each clustering solution compute the dwell time of each PL state
% for all subjects. Perform hypothesis tests to check for differences
% in the dwell time of PL states between conditions.
%
% INPUT:
% data_dir        directory where the LEiDA results are saved
% cond            tags of each condition considered in the experiment
% pair            0 if different subjects in each condition (default); 
%                 1 if same subjects across conditions
% tr              TR of the neuroimaging data
% n_permutations  number of permutations performed in hypothesis tests
% n_bootstraps    number of bootstraps samples within each permutation
%                 sample used in hypothesis tests
%
% OUTPUT:
% LT              dwell time of each FC state in each clustering solution
%                 for all subjects
% LT_pval2sided   two-sided p-value from testing whether the mean dwell time
%                 of a given FC state differs between conditions
% effectsize      Hedge's effect size
% levene_pval     p-value obtained from computing the Levene's test on the
%                 dwell time values
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org



% Default number of permutations
% n_permutations = 10000;
% Default number of bootstrap samples within each permutation sample
% n_bootstraps = 500;

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with K-means results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_V1])
    load([data_dir file_V1], 'Time_sessions', 'Data_info', 'idx_data');
end
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Number of scans considered to compute V1
N_scans = length(idx_data);

% Number of conditions of the experiment
n_Cond = size(cond,2);

Index_Conditions = [];
for s = 1:length(idx_data)
    FileName = Data_info(idx_data(s)).name;
    find_c = 0;
    for c = 1:n_Cond
        if contains(FileName,string(cond(c)))
            find_c = 1;
            Index_Conditions = cat(2, Index_Conditions, c);
        else
            continue
        end
    end
    if find_c == 0
        Index_Conditions = cat(2, Index_Conditions, 0);
    end    
end

% Code below computes the dwell time and runs hypothesis tests
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DWELL TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')
%% CALCULATE THE DWELL TIME OF EACH STATE FOR K-MEANS CLUSTERING SOLUTIONS

LT = zeros(N_scans,length(rangeK),rangeK(end));

for k = 1:length(rangeK)
    for s = 1:N_scans
        
        T = Time_sessions == idx_data(s);
        % Vector of the cluster to which each observation belongs to for subject s
        Ctime = Kmeans_results{k}.IDX(T);
        
        for c = 1:rangeK(k)
            
            Ctime_bin = Ctime == c;
            
            % vector of entries in the FC state under analysis
            in = find(diff(Ctime_bin) == 1);
            % vector of exits from the FC state under analysis
            out = find(diff(Ctime_bin) == -1);
        
            num_periods = 0;
            sum_periods = 0;
            
            % The duration of the FC state active at the beginning and end
            % of each scan is not considered
            if length(out) > length(in)
                out = out(2:end);
            elseif length(in) > length(out)
                in = in(1:end-1);
            elseif ~isempty(in) && ~isempty(out) && out(1) < in(1)
                out = out(2:end);
                in = in(1:end-1);
            end

            if isempty(in) && isempty(out)
                LT(s,k,c) = 0;
            else
                for i=1:length(in)
                    num_periods = num_periods + 1;
                    sum_periods = sum_periods + (out(i) - in(i));
                end
                LT(s,k,c) = (sum_periods / num_periods) * tr;
            end        
        end
    end
end

clear Kmeans_results

%% CHECK HOMOGENEITY OF VARIANCES USING LEVENE'S TEST

if pair == 0
    
    disp('Running Levene''s test of equality of variances to select the statistic to run the permutation tests')
    % matrix to store p-value of Levene's test of equality of variances
    levene_pval = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));

    for k = 1:length(rangeK)
        % disp(['- ' num2str(rangeK(k)) ' FC states'])
        for c = 1:rangeK(k)
            cond_pair = 1;
            for cond1 = 1:n_Cond-1
                for cond2 = cond1+1:n_Cond
                    a = squeeze(LT(Index_Conditions == cond1,k,c));  % Vector containing LT of c in condition 1
                    b = squeeze(LT(Index_Conditions == cond2,k,c));  % Vector containing LT of c in condition 2

                    % data points of dwell time for each clustering
                    % solution from different conditions
                    data_vec = cat(1,a,b);
                    % column vector with index of groups
                    groups = cat(1,zeros(numel(a),1),ones(numel(b),1));
                    % perform Levene's Test of equality of variances (store p-value)
                    levene_pval(cond_pair,k,c) = vartestn(data_vec, groups,'TestType','LeveneAbsolute','Display','off');
                    cond_pair = cond_pair+1;
                end
            end
        end
    end
end
%% PERMUTATION STATISTICS WITH WITHIN BOOTSTRAP SAMPLES

disp(' ')
disp(['Testing intergroup differences in dwell time using ' num2str(n_permutations) ' permutations:'])

% Store p-values for K-means clustering solutions
LT_pval = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));
LT_pval2sided = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));
effectsize = zeros(n_Cond*(n_Cond-1)/2,length(rangeK),rangeK(end));

for k = 1:length(rangeK)
    disp(['- K = ' num2str(rangeK(k))])
    for c = 1:rangeK(k)
        cond_pair = 1;
        for cond1 = 1:n_Cond-1
            for cond2 = cond1+1:n_Cond
                a = squeeze(LT(Index_Conditions == cond1,k,c))';  % Vector containing LT of c in condition 1
                b = squeeze(LT(Index_Conditions == cond2,k,c))';  % Vector containing LT of c in condition 2
                
                if pair == 1
                    % disp([cond{cond1} ' vs ' cond{cond2} ' (using paired-sample t-test statistic)'])
                    stats = bootstrap_within_permutation_paired_samples([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],...,
                                                                        n_permutations,n_bootstraps,0.05);
                else
                    if levene_pval(cond_pair,k,c) < 0.05 % null hypothesis of homogeneity of variances rejected
                        % disp([cond{cond1} ' vs ' cond{cond2} ' (using Welch''s t-test statistic)'])
                        stats = bootstrap_within_permutation_ttest2([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],n_permutations,...,
                                                                     n_bootstraps,0.05,'welchtest');
                    else % null hypothesis of homegeneity of variances not rejected
                        % disp([cond{cond1} ' vs ' cond{cond2} ' (using t-test statistic)'])
                        stats = bootstrap_within_permutation_ttest2([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],n_permutations,...,
                                                                     n_bootstraps,0.05,'ttest');
                    end
                end

                % p-values and effect size
                LT_pval(cond_pair,k,c) = min(stats.pvals);
                LT_pval2sided(cond_pair,k,c) = 2*min(stats.pvals);
                effectsize(cond_pair,k,c) = stats.eff;
                cond_pair = cond_pair + 1;
            end
        end
    end
end

% Name of the file to save output
save_file = 'LEiDA_Stats_DwellTime.mat';

% Save K-means clustering solutions results:
if pair == 0
    save([data_dir '/' save_file],'LT','LT_pval','LT_pval2sided', 'effectsize', 'levene_pval',...,
                              'cond','rangeK','file_cluster','file_V1','Index_Conditions')
else
    save([data_dir '/' save_file],'LT','LT_pval','LT_pval2sided', 'effectsize',...,
                              'cond','rangeK','file_cluster','file_V1','Index_Conditions')
end

disp(' ')                          
disp(['Dwell time values and results from permutation tests saved successfully as ' save_file])
disp(' ')