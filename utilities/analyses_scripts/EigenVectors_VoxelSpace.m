function EigenVectors_VoxelSpace(data_dir,save_dir,leida_dir,selectedK)
%
% This function loads the fMRI signal in each voxel of a brain template.
% Resizes the template to the 10mm MNI space. Computes the fMRI phase
% leading eigenvectors for each TR for all participants. Then computes the
% mean leading eigenvector for each LEiDA state across conditions.
%
% INPUT:
% data_dir      directory where the fMRI data in NIFTI format is saved
% save_dir      directory to save the new data
% leida_dir     directory where the results from running LEiDA are saved
% selectedK     value of K to be further analysed
%
% OUTPUT:
% V1_MNI10mm    leading eigenvectors in MNI 10mm space
% mean_V1       mean leading eigenvectors based on the time courses of PL
%               states obtained from running the K-means algorithm
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% INPUT EXAMPLES:
% data_dir = 'D:/LEiDA_Toolbox/Outputs/dparsf/nofilt_noglobal/func_preproc/';
% save_dir = 'D:/LEiDA_Toolbox/ABIDE_dparsf_MNI10mm/';
% leida_dir = 'D:/LEiDA_Toolbox/LEiDA_Results_ABIDE_dparsf_AAL120/';
% selectedK = 15;

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([leida_dir file_cluster])
    load([leida_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Gather the PL state time courses across all participants
state_time_idxs = Kmeans_results{rangeK == selectedK}.IDX;

% Get number of files in folder
data_info = dir([data_dir '*.nii.gz']);
num_subjs = numel(data_info);

% Load 10mm MNI voxel space
MNI10mm_Mask = niftiread('MNI152_T1_10mm_brain_mask.nii');
sz = size(MNI10mm_Mask); % size of the 10mm MNI mask
ind_voxels = find(MNI10mm_Mask(:) > 0); % find the non-zero elements in the mask
n_voxels = length(ind_voxels);

% Matrix to store the leading eigenvectors of all subjects at each TR
V1_all = zeros(length(state_time_idxs)*2,n_voxels);

t_all = 0;
for s = 1:num_subjs
    
    file = data_info(s).name;
    [~, baseFileName, ~] = fileparts(file);
    ix = strfind(baseFileName,'_'); % get the underscore locations
    if length(ix) == 4 % File with names like CMU_a_0050
        saveFileName = baseFileName(1:ix(3)); % return the substring up to 3rd underscore
    else
        saveFileName = baseFileName(1:ix(2)); % return the substring up to 2nd underscore
    end
    
    if size(file,1)
        disp(['Participant ' saveFileName(1:end-1) ':']);
        
        % Read the nii file
        fMRI_MNI = niftiread([data_dir file]);     
        T = size(fMRI_MNI,4); % number of volumes
        
        disp('- Resizing NIFTI file to MNI 10mm space');
        % Files will be resized in order to be accomodated to the MNI10mm template
        fMRI_MNI10mm = zeros(sz(1), sz(2), sz(3), T);
        for t = 1:T
            fMRI_MNI10mm(:,:,:,t) = imresize3(fMRI_MNI(:,:,:,t),sz);
        end
        clear fMRI_MNI
        
        disp('- Computing the fMRI signal phases using the Hilbert transform');
        % Store the fMRI signal phase using the Hilbert transform
        fMRI_ts = zeros(n_voxels,T);
        
        for v = 1:n_voxels
            [I1,I2,I3] = ind2sub(sz,ind_voxels(v));
            fMRI_ts(v,:) = squeeze(fMRI_MNI10mm(I1,I2,I3,:))';
        end
        clear fMRI_MNI10mm
        
        disp('- Saving fMRI time series in MNI 10mm space');        
        save([save_dir saveFileName 'MNI10mm'], 'fMRI_ts')
        
        disp('- Computing the leading eigenvectors');            
        % De-meaning the fMRI signal
        for v = 1:n_voxels
            fMRI_ts(v,:) = fMRI_ts(v,:) - mean(fMRI_ts(v,:));
        end

        % Get the fMRI signal phase using the Hilbert transform
        for v = 1:n_voxels
            fMRI_ts(v,:) = angle(hilbert(fMRI_ts(v,:)));
        end

        % Compute leading eigenvector of each phase coherence matrix
        for t = 2:T-1 % exclude 1st and last TR of each fMRI signal

            % Save the leading eigenvector for time t
            [v1,~] = eigs(cos(fMRI_ts(:,t) - fMRI_ts(:,t)'),1);

            if sum(v1) > 0 % for eigenvectors with sum of entries > 0
                v1 = -v1;
            end

            t_all = t_all + 1; % time point in V1_all

            % row t_all correponds to the computed eigenvector at time t for subject s
            V1_all(t_all,:) = v1;
        end
    end
end
% Reduce size in case some scans have less TRs than tmax
% In this case, these lines of code will not result in changes
V1_all(t_all+1:end,:) = [];

disp(' ');
disp('- Saving leading eigenvectors in MNI 10mm space');
save([leida_dir  'V1_all_MNI10mm'], 'V1_all', '-v7.3')


% Save the mean PL states in MNI 10mm space
mean_V1 = zeros(selectedK,n_voxels);

disp(' ');
if size(V1_all,1) ~= size(state_time_idxs,2)
    error('- Number of leading eigenvectors and length of state time courses do not coincide');
else
    for k = 1:selectedK
        disp(['- Computing mean leading eigenvector for PL state ' num2str(k)]);
        idx_k = state_time_idxs == k;
        mean_V1(k,:) = mean(V1_all(idx_k,:),1);
    end
end 
disp(' ');
disp(['- Saving the mean leading eigenvectors for K = ' num2str(selectedK)]);
% Create a directory to store results for defined value of K
if ~exist([leida_dir 'K' num2str(selectedK) '/'], 'dir')
    mkdir([leida_dir 'K' num2str(selectedK) '/']);
end
K_dir = [leida_dir 'K' num2str(selectedK) '/'];
save([K_dir  'V1_VoxelSpace'], 'mean_V1', 'ind_voxels', 'MNI10mm_Mask');