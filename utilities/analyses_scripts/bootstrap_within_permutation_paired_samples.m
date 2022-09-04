function [stats] = bootstrap_within_permutation_paired_samples(data,design,niter,nboot,pthr)
% PERMUTATION_HTEST2_NP - A "non parametric" two-sample hypotesis test that, instead of
% relying on the test-type standard distribution, uses permutations of group labels to
% estimate the null distribution. The null distribution is computed
% independently for each data point (= row), i.e. we do not assume the same
% distribution for each datapoint. However, we do assume that the data
% points are comparable (e.g. they correspond to the same location
% collected across all subjects)
%
% INPUT:
% data    a matrix where each column is a subject and each row is a
%         data-point for example a voxel intensity in fMRI, a node level
%         value in a network, etc. NaN values will be ignored.
% design  a row vector containing the numbers 1 and 2 for the two groups
% niter   number of permutations (recommended 5000)
% htest   hypothesis test used to compare populations. The script is
%         prepared to run the ttest2, kstest2, and ranksum tests. 
%
%
% OUTPUT:
% stats is a struct with the following subfields:
% pvals   p-values for each datapoint; it returns in order the p-values
%         for the right tail and for the left tail
% tvals  test statistic values for datapoint, positive tvals mean 
%        group 1 > group 2
%
% Notes: the null distribution is estimated using the matlab function
% ksdensity by interpolating the permuted data. The distribution is
% estimated over 200 points if niter<=5000, otherwise it is estimated over
% round(200*niter/5000) points, for greater precision.

%  Miguel Farinha January 2022 (miguel.farinha@ccabraga.org)
%  adapted from: Enrico Glerean 2013 and Henrique Fernandes 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsubj = size(data,2);   % number of subjects
if(size(design,2) ~= Nsubj)
    error('Mismatched number of subjects: the number of columns of data variable  should match the number of columns of the design variable.')
end
if(size(design,1) ~= 1)
    error('The design variable should only contain 1 row')
end

g1 = find(design==1); % vector with the positions of cond1 in the design vector
g2 = find(design==2); % vector with the positions of cond2 in the design vector
if((length(g1)+length(g2))~=Nsubj)
    error('The design variable should only contain numbers 1 and 2.')
end

if(niter<=0)
    disp('The variable niter should be a positive integer, function will continue assuming niter=5000.')
    niter=5000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYPOTHESIS TESTING (for each row/area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats.tvals=tt_np(data,g1,g2); % similar to ttest2

NC    = size(data,1); % number of comparisons
tvals = zeros(NC,1);
diffs = zeros(NC,1);

for t = 1:NC
    [~,~,~,STATS]  = ttest(data(t,g1)',data(t,g2)',pthr);
    tvals(t,:)     = STATS.tstat;
    diffs(t,:)     = mean(data(t,g1)) - mean(data(t,g2));
end
    
stats.tvals = tvals;
% tvals(isnan(tvals)) = 0;   % or tvals(tvals~=tvals) = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERMUTATION TESTING (for each row/area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outptus the pval (from the computed null distribution using permutation 
% testing) given the tstat previously calculated.
% each comparison is treated independently

%changes this to work for median values close to 0 
pvals = zeros(NC,2);

parfor n=1:NC
    if sum(abs(data(n,g1)))>0 || sum(abs(data(n,g2)))>0  % Exclude tests where all (tstat=NaN) or most of the population (median=0) as a null value.
        pvals(n,:) = tt_np_pval(data(n,:),g1,g2,niter,nboot,tvals(n));
     else
        disp('Data with lots of zeros!');
        pvals(n,:) = [NaN NaN];
     end
end

stats.pvals = pvals;

% two-sided p-value:
stats.pvals_2sided = 2*min(pvals);

% Hedge's effect size:
n1 = length(g1);
n2 = length(g2);
m1 = mean(data(g1));
m2 = mean(data(g2));
s1 = std(data(g1));
s2 = std(data(g2));
s = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
eff = abs(m1-m2)/s;

stats.eff = eff;

stats.diffs = diffs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pval = tt_np_pval(data,g1,g2,niter,nboot,tval)
    
    outiter = zeros(niter,1);
    ND = length(data);
    for iter = 1:niter
        perm = randperm(ND); % random permutation of the indices of the subjects' ID
        temp = data(perm); % temp is the vector of data with positions changed according to perm
        
        % (1) get sample sizes (size of both conditions is equal)
        xnans = isnan(temp(:,g1));
        if any(xnans(:))
            nx = sum(~xnans,2);
        else
            nx = size(temp(:,g1),2); 
        end
        
        % (2) compute the mean of the difference of the permuted samples
        temp_g1 = temp(:,g1);
        temp_g2 = temp(:,g2);
        d = temp_g1 - temp_g2;
        mean_d = nanmean(d);
        
        % (3) get nboot bootstrap samples from permuted sample to compute
        % variance of each group
        [~,bootsam_g1] = bootstrp(nboot, [], g1);
        [~,bootsam_g2] = bootstrp(nboot, [], g2);
        
        % (4) vectors to store the standard error for each bootstrap sample
        se_boot = zeros(nboot,1);
        for b = 1:nboot
            % get bootstrap sample for each group
            tempboot_g1 = temp_g1(bootsam_g1(:,b)); % dim = 1*len(g1)
            tempboot_g2 = temp_g2(bootsam_g2(:,b)); % dim = 1*len(g2)
            % compute the bootstrap sample differences
            d_boot = tempboot_g1 - tempboot_g2;
            % compute the mean of bootstrap differences
            dm_boot = nanmean(d_boot);
            % compute the standard error of differences
            se_boot(b) = sqrt((sum((d_boot-dm_boot).^2))/(nx*(nx-1)));
        end
        
        % (5) compute the mean of the standard error to be used to compute
        % the value of the statistic under the current permutation sample
        se_perm = nanmean(se_boot);
        
        % (6) compute the value of the statistic
        ratio = mean_d / se_perm;
        if ~isnan(ratio)
            outiter(iter) = ratio;
        end
    end
    
    NCDF = 200;
    if(niter > 5000)
        NCDF = round(200*niter/5000);
    end
    [fi, xi] = ksdensity(outiter,'function','cdf','npoints',NCDF); % estimated cumulative distribution function

    % trick to avoid NaNs, we approximate the domain of the CDF between
    % -Inf and Inf using the atanh function and the eps matlab precision
    % variable
    
    pval_left = interp1([atanh(-1+eps) xi atanh(1-eps)],[0 fi 1],tval); % G1 > G2
    pval_right = 1-pval_left; % G1 < G2
    pval = [pval_right pval_left];
end
 