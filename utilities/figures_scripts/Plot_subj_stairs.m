function Plot_subj_stairs(data_dir,save_dir,selectedK,subj)
%
% Plot the state time course for a specific participant considering K PL
% states as a stairs plot
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% subj          subject defined by the user
%
% OUTPUT:
% .fig/.png     Plot of the state time course for a given value of K and
%               subject as a stairs plot
%
% Authors: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_V1])
    load([data_dir file_V1], 'Time_sessions', 'Data_info', 'idx_data');
end
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Get the fMRI signal from the subject of interest
for i = 1:length(idx_data)
    s = idx_data(i);
    if contains(Data_info(s).name,string(subj))
        signal = importdata([Data_info(s).folder '/' Data_info(s).name]);
        subj_name = Data_info(s).name;
        if isstruct(signal)
        % fMRI data should be stored in a field called data
            signal = signal.data;
        end
        T = Time_sessions == s;
    end
end

% Get Tmax from the signal
[~,Tmax] = size(signal);

% Getting cluster assignments for the specific subject
Ctime = Kmeans_results{rangeK == selectedK}.IDX(T);

% Colormap for the figure
if selectedK < 13
    cmap = linspecer(selectedK, 'qualitative');
else
    cmap = linspecer(selectedK);
end

X = linspace(2,Tmax-1,Tmax-2)';
Y = zeros(Tmax-2,selectedK+1);
for i = 1:size(Y,2)-1
    Y(:,i) = (selectedK+1-i)+0.5*(Ctime == i);
end

disp(' ');
disp(['Plotting stairs state time courses for ' num2str(selectedK) ' clusters for participant ' subj_name ':'])
Fig = figure('Position', get(0, 'Screensize'));
subplot(1,1,1)
h = stairs(X,Y,'LineWidth',1.5);
for j = 1:selectedK
    h(j).Color = cmap(j,:);
end
h(selectedK + 1).Color = 'black';
h(selectedK + 1).LineStyle = 'none';
legend_titles = cell(1,selectedK+1);
legend_titles{1} = '';
for j = 1:selectedK
    legend_titles{j+1} = ['PL state ' num2str(selectedK +1 - j)];
end
yticklabels(legend_titles);
xlim([1 Tmax])
xticks(2:20:Tmax-2);
xlabel('Time (TR)', 'Fontsize', 12);
set(gca,'TickLength',[0 0])
box off
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterStairs.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterStairs.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_ClusterStairs']);

close all;
