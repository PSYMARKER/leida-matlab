function [V1_all,Time_sessions,Data_info] = LEiDA_data(data_dir,save_dir,n_areas,tmax,filter,flp,fhi,tr)
%
% For each subject compute the leading eigenvetor of the phase coherence
% matrix calculated at each recording frame.
%
% INPUT:
% data_dir      directory where the parcellated fMRI data are saved;
%               fMRI data should be one file per subject; data can be a
%               matrix containing the time series in the formats .mat,
%               .1D and .txt; if data is as a struct then fMRI signal
%               should correspond to a field called data
% save_dir      directory where the leading eigenvectors will be saved
% n_areas       number of brain areas to consider for analysis
% tmax          maximum number of volumes across fMRI sessions
% filter        0, temporal filtering (default); 1, no temporal filtering
% flp           lowpass frequency of filter
% fhi           highpass frequency of filter
% tr            TR of fMRI data
%
% OUTPUT:
% V1_all        (n_scans*(tmax-2) x n_areas) leading eigenvectors of all
%                subjects at each time point
% Time_sessions (1 x n_scans*(tmax-2)) scan number of each leading
%                eigenvector
% Data_info      parcellated data
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEADING EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

% Get the files containg the data (.mat or .1D or .txt)
Data_info = [dir([data_dir '*.mat']); dir([data_dir '*.1D']); dir([data_dir '*.txt'])];
% Total number of scans that will be read
n_scans = size(Data_info,1);
disp(['Total number of scans in folder: ' num2str(n_scans)])
disp(' ')

% Matrix to store the leading eigenvectors of all subjects at each TR
V1_all = zeros(n_scans*(tmax-2),n_areas);
% Row vector with the scan number of each leading eigenvector
Time_sessions = zeros(1,n_scans*(tmax-2));  

t_all = 0;
discarded = [];
idx_data = [];
for s = 1:n_scans
    
    disp(['Computing the leading eigenvectors for scan ' num2str(s) ' ' Data_info(s).name]);
    
    % Handling the data differently depending on the file extension:
    % Regardless of the extension we import the data
    signal = importdata([data_dir Data_info(s).name]);
    if isstruct(signal)
        % Load the fMRI signal from the struct signal variable
        % fMRI data should be stored in a field called data
        signal = signal.data;
    end
    % Selecting only the areas specified by the user
    signal = signal(1:n_areas,:);
    
    if any(isnan(signal(:))) || any(isinf(signal(:))) || any(all(signal == 0,2))
        disp('      - NaN, Inf or rows of 0s were found -> Discarded');
        discarded = cat(2,discarded,s);
    else
        % Add index of participant to idx_data
        idx_data = cat(2,idx_data,s);
        
        % De-meaning the fMRI signal
        for n = 1:n_areas
            signal(n,:) = detrend(signal(n,:) - mean(signal(n,:)));
        end
        
        % Apply temporal filtering to the parcellated data
        if filter
            signal = TemporalFiltering(signal,flp,fhi,tr);
        end

        % Get the fMRI signal phase using the Hilbert transform
        for seed = 1:n_areas
            signal(seed,:) = angle(hilbert(signal(seed,:)));
        end

        % Compute leading eigenvector of each phase coherence matrix
        for t = 2:size(signal,2)-1 % exclude 1st and last TR of each fMRI signal

            % Save the leading eigenvector for time t
            [v1,~] = eigs(cos(signal(:,t)-signal(:,t)'),1);

            if sum(v1) > 0 % for eigenvectors with sum of entries > 0
                v1 = -v1;
            end

            t_all = t_all + 1; % time point in V1_all

            % row t_all correponds to the computed eigenvector at time t for subject s
            V1_all(t_all,:) = v1;

            % to which subject the computed leading eigenvector belongs
            Time_sessions(t_all) = s;
        end
    end
end
% Reduce size in case some scans have less TRs than tmax
% In this case, these lines of code will not result in changes
V1_all(t_all+1:end,:) = [];
Time_sessions(:,t_all+1:end) = [];

disp(' ')
disp(['Total number of scans used to compute the leading eigenvectors: ' num2str(length(idx_data))]);

for d = 1:length(discarded)
    disp(['Attention: ' Data_info(discarded(d)).name ' was discarded'])
end
disp(' ')

% Name of the file to save output
save_file = 'LEiDA_EigenVectors.mat';

save([save_dir save_file], 'V1_all','Time_sessions','Data_info','idx_data')

disp(['fMRI Phase Leading Eigenvectors saved successfully as ' save_file])
disp(' ')