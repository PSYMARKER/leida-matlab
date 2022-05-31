function Parcellate_ABIDE_func_preproc(data_dir,info_file,save_dir,parcellation,n_areas)
%
% Parcellate the ABIDE functional preprocessed files and save the files
% with the tag of the condition appended to the file name.
%
% INPUT:
% data_dir      directory where the .nii.gz files with the ABIDE data are
%               stored
% info_file     path to file with the phenotypic data from each subject
% save_dir      directory to save the new data
% parcellation  parcellation atlas to be used to parcel the data
% n_areas       number of areas of the parcellation to consider
%
% Author: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%         Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% INPUT EXAMPLES:
% data_dir = 'D:/LEiDA_Toolbox/Outputs/dparsf/nofilt_noglobal/';
% info_file = 'D:/LEiDA_Toolbox/Phenotypic_V1_0b.csv';
% save_dir = 'D:/LEiDA_Toolbox/ABIDE_dparsf_AAL120/';
% parcellation = 'AAL120'; (AAL116, AAL120, dbs80, glasser360)
% n_areas = 94; (excluding cerebellum)

% Get number of files in folder
aux_data = dir([data_dir '*.nii.gz']);
num_subjs = numel(aux_data);

% Order the directory by name
[~,ind] = sort({aux_data.name});
data_info = aux_data(ind);

% Read the CSV with phenotypic data
info_data = readtable(info_file);
len_info = size(info_data,1);

% Obtain parcellation atlas from ParcelsMNI2mm
V_Parcel = struct2array(load('ParcelsMNI2mm',['V_' parcellation])); 
sz = size(V_Parcel); % size of the parcellation

% Count the number of subjects and maximum number of TRs
n_hc = 0;
n_ad = 0;
n_asp = 0;
n_pdd = 0;
tmax = 0;
for s = 1:num_subjs
    
    file = data_info(s).name;
    [~, baseFileName, ~] = fileparts(file);
    
    if size(file,1)
        disp(['Parcellating data file ' baseFileName ' using parcellation ' parcellation]);
        
        % Read the nii file
        fMRI_MNI = niftiread([data_dir file]);     
        T = size(fMRI_MNI,4); % number of volumes
        
        % Check if nii files correspond to MNI2mm
        if size(fMRI_MNI,1) ~= sz(1) || size(fMRI_MNI,2) ~= sz(2) || size(fMRI_MNI,3) ~= sz(3)
            % Files will be resized in order to be accomodated to the MNI2mm template
            fMRI_MNI2mm = zeros(sz(1), sz(2), sz(3), T);
            for t = 1:T
                fMRI_MNI2mm(:,:,:,t) = imresize3(fMRI_MNI(:,:,:,t),sz);
            end
        end
                
        data = zeros(n_areas,T);
        
        for n = 1:n_areas
            ind_voxels = find(V_Parcel == n);
            
            for v = 1:numel(ind_voxels)
                [I1,I2,I3] = ind2sub(sz,ind_voxels(v));
                if ~isnan(fMRI_MNI2mm(I1,I2,I3,1))
                    data(n,:) = data(n,:) + squeeze(fMRI_MNI2mm(I1,I2,I3,:))';
                end
            end

            data(n,:) = data(n,:)/numel(ind_voxels);
            data(n,:) = detrend(data(n,:) - mean(data(n,:)));
            data(n,:) = data(n,:)/std(data(n,:));
        end
        
        for i = 1:len_info
            pattern = num2str(info_data{i,2});
            if contains(baseFileName,pattern)
                
                ix = strfind(baseFileName,'_'); % get the underscore locations
                if length(ix) == 4 % File with names like CMU_a_0050
                    saveFileName = baseFileName(1:ix(3)); % return the substring up to 3rd underscore
                else
                    saveFileName = baseFileName(1:ix(2)); % return the substring up to 2nd underscore
                end
                
                if tmax < size(data,2)
                    tmax = size(data,2);
                end
                
                % If DSM_IV_TR not specified check if DX_GROUP == Control
                if isequal(info_data{i,4},-9999) && isequal(info_data{i,3},2)
                    n_hc = n_hc + 1;
                    disp('          - Participant -> CONTROL (DX_GROUP)');
                    save([save_dir saveFileName parcellation '_CONT'], 'data')
                % If DSM_IV_TR not specified check if DX_GROUP == Autism
                elseif isequal(info_data{i,4},-9999) && isequal(info_data{i,3},1)
                    n_ad = n_ad + 1;
                    disp('          - Participant -> AUTISM (DX_GROUP)');
                    save([save_dir saveFileName parcellation '_AUT'], 'data')
                % DSM_IV_TR: Control == 0    
                elseif isequal(info_data{i,4},0)
                    n_hc = n_hc + 1;
                    disp('          - Participant -> CONTROL');
                    save([save_dir saveFileName parcellation '_CONT'], 'data')
                % DSM_IV_TR: Autism == 1
                elseif isequal(info_data{i,4},1)
                    n_ad = n_ad + 1;
                    disp('          - Participant -> AUTISM');
                    save([save_dir saveFileName parcellation '_AUT'], 'data')
                % DSM_IV_TR: Aspergers == 2
                elseif isequal(info_data{i,4},2)
                    n_asp = n_asp + 1;
                    disp('          - Participant -> ASPERGERS');
                    save([save_dir saveFileName parcellation '_ASP'], 'data')
                % DSM_IV_TR: PDD-NOS == 3
                elseif isequal(info_data{i,4},3)
                    n_pdd = n_pdd + 1;
                    disp('          - Participant -> PDD-NOS');
                    save([save_dir saveFileName parcellation '_PDD_NOS'], 'data')
                % DSM_IV_TR: Aspergers or PDD-NOS == 4 (only 6 subjects)
                elseif isequal(info_data{i,4},4)
                    n_ad = n_ad + 1;
                    disp('          - Participant -> ASPERGERS OR PDD-NOS');
                    save([save_dir saveFileName parcellation '_AUT'], 'data')
                end
            end
        end
    end
end
disp(['Number of participants with tag Control: ' num2str(n_hc)]); 
disp(['Number of participants with tag Autism: ' num2str(n_ad)]); 
disp(['Number of participants with tag Aspergers: ' num2str(n_asp)]); 
disp(['Number of participants with tag PDD-NOS: ' num2str(n_pdd)]); 
disp(['The maximum number of TRs across participants is: ' num2str(tmax)]);

