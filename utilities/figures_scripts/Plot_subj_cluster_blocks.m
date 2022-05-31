function Plot_subj_cluster_blocks(data_dir,save_dir,selectedK,subj)
%
% Plot the state time course for a specific participant considering K PL
% states along with the fMRI signal
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% subj          subject defined by the user
%
% OUTPUT:
% .fig/.png     Plot of the state time course for a given value of K and
%               subject along with the fMRI signal
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

% Get the number of areas considered for analysis
n_areas = size(Kmeans_results{rangeK == selectedK}.C,2);

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
        signal = signal(1:n_areas,:);
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

disp(' ');
disp(['Plotting state cluster blocks for K =' num2str(selectedK) ' clusters for participant ' subj_name ':'])
Fig = figure('Position', get(0, 'Screensize'));
subplot(2,1,2)
ymax = max(abs(signal(:)));
hold on
for t = 2:Tmax-1
    x = [t t+1 t+1 t];
    x = (x-1); %.*TR
    y = [-ymax -ymax ymax ymax];
    p = patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',cmap(Ctime(t-1),:),'FaceAlpha',0.9);
end
plot(2:1:Tmax-1,signal(:,2:Tmax-1),'color',[.5 .5 .5])
hold on
patches = gobjects(1,selectedK);
for i = 1:selectedK
    patches(i) = fill(nan, nan, cmap(i,:),'FaceAlpha',0.9);
end
legend_titles = cell(1,selectedK);
for j = 1:selectedK
    legend_titles{j} = ['PL state ' num2str(j)];
end
legend(patches,legend_titles,'Location','northoutside','Orientation','horizontal','Box','off','Fontsize',10);
ylabel('fMRI signal','Fontsize',10);
xlabel('Time (TR)','Fontsize',10)
xlim([2 Tmax-1])
ylim([-ymax ymax])
xticks(2:25:Tmax-1)

saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterBlocks.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_ClusterBlocks.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_ClusterBlocks']);

close all;



