function Plot_K_state_time(data_dir,save_dir,selectedK)
%
% Plot the state time courses for all participants and make pie plots with
% the percentage of time spent in each state.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
%
% OUTPUT:
% .fig/.png     Plot of the state time courses for all subjects by
%               condition with pie plots of fractional occupancy of states
%
% Authors: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

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
    load([data_dir file_P], 'cond', 'Index_Conditions');
end

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Compute the Tmax from the Time_sessions variable
[Ts,~] = hist(Time_sessions,unique(Time_sessions));
Tmax = max(Ts(:));

% Compute the number of subjects belonging to each condition
[nSub_cond,~] = hist(Index_Conditions,unique(Index_Conditions));
nSub_cond_max = max(nSub_cond(:));

% Vector to save the state time courses for each subject by condition
State_time_cond = zeros(nSub_cond_max,Tmax,n_Cond);
for c = 1:n_Cond
    subs_cond = find(Index_Conditions == c);
    for s = 1:length(subs_cond)
        subj = subs_cond(s);
        T = Time_sessions == idx_data(subj);
        Tmax_s = length(find(T == 1));
        State_time_cond(s,1:Tmax_s,c) = Kmeans_results{rangeK == selectedK}.IDX(T);
    end
end
% Replace 0s by NaN (subjects may have less TRs than Tmax):
State_time_cond(State_time_cond == 0) = NaN;

% Select most appropriate sliding-window to compute states' fractional occupancy
slide_windows = cell(1,7);
n_windows = [4,5,6,7,8,9,10];
last_point = zeros(1,7);
for i = 1:length(n_windows)
    window = floor(size(State_time_cond,2)/n_windows(i));
    slide_windows{i} = 1:window:size(State_time_cond,2);
    last_point(i) = slide_windows{i}(end);
end
[~,ind_max] = max(last_point);
slide_points = slide_windows{ind_max};
slide_points(end) = size(State_time_cond,2);

disp(' ');
disp(['Plotting state time coures obtained for ' num2str(selectedK) ' clusters for all participants by condition:'])
Fig = figure('Position', get(0, 'Screensize'));
for c = 1:n_Cond
    subplot(n_Cond*2,length(slide_points),[1+(c-1)*length(slide_points)*2:(c*2-1)*length(slide_points)])
    colormap(linspecer)
    h = imagesc(State_time_cond(1:nSub_cond(c),:,c));
    set(h, 'AlphaData', ~isnan(State_time_cond(1:nSub_cond(c),:,c))) % plot NaNs with background colour
    set(gca,'TickLength',[0 0])
    Ylm = ylim; % get x, y axis limits 
    Xlm = xlim;                          % so can position relative instead of absolute
    Xlb = 0.99*max(Xlm);                    % set horizontally at midpoint
    Ylb = 0.99*max(Ylm);                  % and just 1% below minimum y value
    hXLbl = xlabel('Time','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 
    ylabel('Subjects')
    title(cond{c},'interpreter','none')
    box off
    hold on
    plot([slide_points(2:end-1); slide_points(2:end-1)],[0*ones(1,length(slide_points)-2); 2*nSub_cond(c)*ones(1,length(slide_points)-2)],'k','linewidth',1.5)
    
    for slide = 1:length(slide_points)-1
        State_Cond = State_time_cond(:,slide_points(slide):slide_points(slide)+window-1,c);
        subplot(n_Cond*2,length(slide_points)-1,(c*2-1)*(length(slide_points)-1)+slide)
        p = pie(histcounts(State_Cond(~isnan(State_Cond)),1:selectedK+1));
        for k = 1:selectedK
            t = p(k*2);
            t.FontSize = 6;
        end
    end
end
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_StateTimes.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_StateTimes.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_StateTimes']);

close all;
    

