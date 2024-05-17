function LEiDA_stats_TransitionMatrix(data_dir,save_dir,selectedK,n_permutations,n_bootstraps)
%
% For the selected K compute the transition probability matrix for each
% participant. Perform hypothesis tests to check for differences in the
% state-to-state transition probabilities between conditions.
%
% INPUT:
% data_dir        directory where the LEiDA results are saved
% save_dir        directory to save the results from the hypothesis tests
%                 of the comparison of the mean transition probabilities for
%                 the selected K
% selectedK       K defined by the user
% n_permutations  number of permutations performed in hypothesis tests
% n_bootstraps    number of bootstraps samples within each permutation
%                 sample used in hypothesis tests
%
% OUTPUT:
% TM              transition matrix for optimal K for all participants
% TMnorm          transition probability matrix (TPM) for optimal K for all
%                 participants
% TM_pval2sided   two-sided p-value from testing whether the mean transition
%                 probability of a given FC state differs between conditions      
% effectsize      Hedge's effect size
% levene_pval     p-value obtained from computing the Levene's test on the
%                 state-to-state transition probability values
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org


% Default number of permutations
% n_permutations = 10000;
% Default number of bootstrap samples within each permutation sample
% n_bootstraps = 500;

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
if isfile([data_dir file_V1])
    load([data_dir file_V1], 'Time_sessions', 'idx_data');
end
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end
if isfile([data_dir file_P])
    load([data_dir file_P], 'cond', 'P', 'Index_Conditions', 'pair');
end

% Number of scans considered to compute V1
N_scans = length(idx_data);

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Code below computes the fractional occupancy and runs hypothesis tests
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITION PROBABILITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

%% CALCULATE THE TRANSITION MATRIX FOR THE OPTIMAL K

% Matrix with probability of occurrence of each FC state for chosen K
P_K = squeeze(P(:,selectedK-1,1:selectedK));

% Transition probability matrix
TM = zeros(N_scans, selectedK, selectedK);
TMnorm = zeros(N_scans, selectedK, selectedK);

for s = 1:N_scans

    T = Time_sessions == idx_data(s);
    Ctime = Kmeans_results{rangeK == selectedK}.IDX(T);
    % alpha corresponds to departure state;
    for alpha = 1:selectedK
        % beta corresponds to arrival state;
        for beta = 1:selectedK
            % count number transition from alpha to beta
            alpha2beta = 0;
            for i = 1:length(Ctime)-1
                % if alpha = beta then we have 1 transition
                if isequal(Ctime(i),alpha) && isequal(Ctime(i+1),beta)
                    alpha2beta = alpha2beta + 1;
                end
            end
            TM(s,alpha,beta) = alpha2beta;
        end
    end
    % normalise by T-1 as in Vohryzek et al., (2020)
    TM(s,:,:) = TM(s,:,:)./(size(Ctime,2)-1);
    TMnorm(s,:,:) = squeeze(TM(s,:,:))./P_K(s,:)';
end

% Set all Nan values to zero
% If a FC state has probability of occurence 0 then there should be a
% probability of 0 to transition to that state
% numNaN = sum(isnan(TMnorm),'all');
% TMnorm(isnan(TMnorm)) = 0;

clear Kmeans_results

%% CHECK HOMOGENEITY OF VARIANCES USING LEVENE'S TEST

if pair == 0
    
    disp('Running Levene''s test of equality of variances to select the statistic to be used in the permutation tests')
    % matrix to store p-value of Levene's test of equality of variances
    levene_pval = zeros(n_Cond*(n_Cond-1)/2,selectedK,selectedK);

    for c_out = 1:selectedK
        for c_in = 1:selectedK
            % disp(['- Transition from FC state ' num2str(c_out) ' to FC state ' num2str(c_in)])
            cond_pair = 1;
            for cond1 = 1:n_Cond-1
                for cond2 = cond1+1:n_Cond
                    a = squeeze(TMnorm(Index_Conditions == cond1,c_out,c_in));
                    b = squeeze(TMnorm(Index_Conditions == cond2,c_out,c_in));
                    
                    data_vec = cat(1,a,b);
                    % column vector with index of groups, SZ = 0 and HC = 1
                    groups = cat(1,zeros(numel(a),1),ones(numel(b),1));
                    % perform Levene Test of equality of variances (store p-value)
                    levene_pval(cond_pair,c_out,c_in) = vartestn(data_vec,groups,'TestType','LeveneAbsolute','Display','off');
                    cond_pair = cond_pair + 1;
                end
            end
        end
    end
end

%%  PERMUTATION STATISTICS WITH WITHIN BOOTSTRAP SAMPLES

disp(' ')
disp(['Testing intergroup differences in transition probabilities using ' num2str(n_permutations) ' permutations:'])

% Store p-values and effect size
TM_pval = zeros(n_Cond*(n_Cond-1)/2,selectedK,selectedK);
TM_pval2sided = zeros(n_Cond*(n_Cond-1)/2,selectedK,selectedK);
effectsize = zeros(n_Cond*(n_Cond-1)/2,selectedK,selectedK);

for c_out = 1:selectedK 
    for c_in = 1:selectedK
        disp(['- From PL state ' num2str(c_out) ' to PL state ' num2str(c_in)])
        cond_pair = 1;
        for cond1 = 1:n_Cond-1
            for cond2 = cond1+1:n_Cond
                a = squeeze(TMnorm(Index_Conditions == cond1,c_out,c_in))';
                a = rmmissing(a);
                b = squeeze(TMnorm(Index_Conditions == cond2,c_out,c_in))';  % Vector containing Prob of c in condition 2
                b = rmmissing(b);

                if pair == 1
                    % disp([cond{cond1} ' vs ' cond{cond2} ' (using paired-sample t-test statistic)'])
                    stats = bootstrap_within_permutation_paired_samples([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],...,
                                                                        n_permutations,n_bootstraps,0.05);
                else      
                    if levene_pval(cond_pair,c_out,c_in) < 0.05 % null hypothesis of homogeneity of variances rejected
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
                TM_pval(cond_pair,c_out,c_in) = min(stats.pvals);
                TM_pval2sided(cond_pair,c_out,c_in) = 2*min(stats.pvals);
                effectsize(cond_pair,c_out,c_in) = stats.eff;
                cond_pair = cond_pair + 1;
            end
        end
    end
end

% Name of the file to save output
save_file = 'LEiDA_Stats_TransitionMatrix.mat';

% Save K-means clustering solutions results:
if pair == 0
    save([save_dir '/' save_file],'TM','TMnorm','TM_pval','effectsize','levene_pval', 'TM_pval2sided',...,
                              'cond','rangeK','file_cluster','file_V1','Index_Conditions')
else
    save([save_dir '/' save_file],'TM','TMnorm','TM_pval','effectsize','TM_pval2sided',...,
                              'cond','rangeK','file_cluster','file_V1','Index_Conditions')
end

disp(' ')
disp(['Transition probability values and results from permutation tests saved successfully as ' save_file])
disp(' ')
