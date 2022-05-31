function append_tag(data_dir,info_file,tag_col,save_dir)
%
% Example script to append the condition tag of each participant to the
% filename to % enable running the script LEiDA_Start.
%
% INPUT:
% data_dir      directory where data files are stored
% info_file     path to file with the phenotypic data from each subject
% tag_col       number of the column which determines the condition of
%               each participant
% save_dir      directory to save the new data
%
% Author: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% Input examples:
% data_dir = 'D:\LEiDA_Toolbox\Outputs\niak\nofilt_noglobal\rois_aal\';
% info_file = 'D:\LEiDA_Toolbox\Phenotypic_V1_0b.csv';
% tag_col = 3;
% save_dir = 'D:\LEiDA_Toolbox\ABIDE_niak_aal_data\';

% Get the files containg the data (.mat or .1D or .txt)
dataFiles = [dir(fullfile([data_dir '*.mat'])); dir(fullfile([data_dir '*.1D'])); dir(fullfile([data_dir '*.txt']))];
fileNames = {dataFiles.name};

% Read the CSV with phenotypic data
info_data = readtable(info_file);
% Get the unique values in the column tag_col
C = unique(info_data{:,tag_col});

for k = 1:length(inputFiles)
    
    [~,thisFileName,ext] = fileparts(fileNames{k});
    % Prepare the input filename.
    inputFullFileName = fullfile(data_dir, [thisFileName ext]);
    % Prepare the output filename.
    for i = 1:length(C)
        if isequal(info_data{s,tag_col},C(i))
            if isa(C,'double')
                outputBaseFileName = sprintf('%s_%s%s',thisFileName,num2str(C(i)),ext);
            elseif isa(C,'cell')
                outputBaseFileName = sprintf('%s_%s%s',thisFileName,C(i),ext);
            end
        end
    end
    outputFullFileName = fullfile(save_dir, outputBaseFileName);
    % Do the copying and renaming all at once.
    copyfile(inputFullFileName, outputFullFileName);
end
