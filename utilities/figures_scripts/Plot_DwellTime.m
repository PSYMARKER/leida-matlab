function Plot_DwellTime(data_dir)
%
% Plot the results from the hypothesis tests obtained from comparing the
% mean dwell time between conditions
%
% INPUT:
% data_dir      directory where the results from running the hypothesis
%               tests on the dwell time of PL states
%
% OUTPUT:
% Fig1          plot of the two-sided p-values obtained across K from the
%               comparison of the mean dwell time of PL states for each
%               pair of conditions
% Fig2          barplot of the mean dwell time of each PL state for all
%               conditions
% Fig3          barplot of the mean dwell time of each PL state for each
%               pair of conditions
% Fig4          plot of the Hedge's effect size obtained across K from the
%               comparison of the mean dwell time of PL states for each
%               pair of conditions
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with results for the dwell time (output from LEiDA_stats_DwellTime.m)
file_LT = 'LEiDA_Stats_DwellTime.mat';

% Load required data:
load([data_dir file_LT],'cond', 'LT', 'LT_pval2sided', 'Index_Conditions', 'effectsize', 'rangeK');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Create a map to place the figures
pos = power(n_Cond-1,2);
index_fig = reshape(1:pos, n_Cond-1, n_Cond-1).';
subplot_map = ones(n_Cond-1);
subplot_map = triu(subplot_map).';

% Define indices for conditions in subplot
subplot_indices = find(subplot_map);

% Code below generates the figures relative to the dwell time
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DWELL TIME FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Possible pairs of conditions comparisons
condRow = zeros(1,n_Cond*(n_Cond-1)/2);
condCol = zeros(1,n_Cond*(n_Cond-1)/2);
cond_pair = 1;
disp('Possible pairs of condition comparisons:')
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        disp([num2str(cond{cond1}) ' - ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

% Plot of two-sided p-values from dwell time hypothesis tests
Fig1 = figure('Position', get(0, 'Screensize'));
disp(' ')
disp('Plotting two-sided p-values from hypothesis (permutation) tests:')
for s_ind = 1:length(subplot_indices)
    
    % disp(['- ' cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}])
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1.5)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1.5)
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1.5)
    
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            if LT_pval2sided(s_ind,k,c) > 0.05
                semilogy(rangeK(k),LT_pval2sided(s_ind,k,c),'.k','Markersize',4);
                % text(rangeK(k),LT_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if LT_pval2sided(s_ind,k,c) < 0.05 && LT_pval2sided(s_ind,k,c) > (0.05/rangeK(k))
                semilogy(rangeK(k),LT_pval2sided(s_ind,k,c),'+r','Markersize',4);
                % text(rangeK(k),LT_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if LT_pval2sided(s_ind,k,c) < (0.05/rangeK(k)) && LT_pval2sided(s_ind,k,c) > (0.05/sum(rangeK))
                semilogy(rangeK(k),LT_pval2sided(s_ind,k,c),'og','Markersize',4);
                % text(rangeK(k),LT_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
            if LT_pval2sided(s_ind,k,c) <= (0.05/sum(rangeK))
                semilogy(rangeK(k),LT_pval2sided(s_ind,k,c),'*b','Markersize',4);
                text(rangeK(k),LT_pval2sided(s_ind,k,c),[' ' num2str(c)])
            end
        end
    end
     
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel([{'Two-sided {\itp}-value'},{'Dwell Time'}])
    xlabel('Number of PL States K')
    xticks([2 4 6 8 10 12 14 16 18 20])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    ylim([10.^(-7.5) 1])
    set(gca,'YTick',10.^(-7:2:0))
end
saveas(Fig1, fullfile(data_dir, 'DwellTime_pvalues.png'),'png');
saveas(Fig1, fullfile(data_dir, 'DwellTime_pvalues.fig'),'fig');
disp('- Plot successfully saved as DwellTime_pvalues');

% Barplot of mean dwell time for each condition across K (all)
if n_Cond > 2
    Fig2 = figure('Position', get(0, 'Screensize'));
    disp(' ')
    disp('Plotting barplot of the mean dwell time across conditions and K:')
    for k = 1:length(rangeK)

        % disp(['- ' num2str(rangeK(k)) ' PL states'])

        for c = 1:rangeK(k)

            subplot_tight(length(rangeK),rangeK(end),c+(k-1)*rangeK(end),0.015)

            hold on

            LT_cond = cell(1,n_Cond);
            mean_LT_cond = zeros(1,n_Cond);
            ste = zeros(1,n_Cond);
            for j = 1:n_Cond
                LT_cond{j} = LT(Index_Conditions == j,rangeK == k+1,c);
                mean_LT_cond(j) = nanmean(LT(Index_Conditions == j,rangeK == k+1,c));
                ste(j) = std(LT(Index_Conditions == j,rangeK == k+1,c))/sqrt(numel(LT(Index_Conditions == j,rangeK == k+1,c)));
            end

            bar(1:n_Cond,mean_LT_cond,'EdgeColor','k','LineWidth',0.8,'FaceColor','none')
            hold on
            errorbar(mean_LT_cond,ste,'LineStyle','none','Color','k','CapSize',5,'LineWidth',0.4);
            if k == length(rangeK)
                set(gca,'XTick',1:n_Cond,'XTickLabel',cond,'Fontsize',6,'TickLabelInterpreter','none')
                % xtickangle(45)
            else
                set(gca,'XTick',1:n_Cond,'XTickLabel',[],'Fontsize',6,'TickLabelInterpreter','none')
            end
            set(gca,'color','none')

            hold off
            box off
        end
    end
    saveas(Fig2, fullfile(data_dir, 'DwellTime_Barplot_Allconditions.png'),'png');
    saveas(Fig2, fullfile(data_dir, 'DwellTime_Barplot_Allconditions.fig'),'fig');
    disp('- Plot successfully saved as DwellTime_Barplot_Allconditions');
end

% Barplot of mean dwell time for each pair of conditions
disp(' ')
disp('Plotting barplots of mean dwell time for each pair of conditions:')
n_compare = n_Cond*(n_Cond-1)/2;
for i = 1:n_compare
    Fig3 = figure('Position', get(0, 'Screensize'));
    % disp(['- ' cond{condRow(i)} ' vs ' cond{condCol(i)}])
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            
            subplot_tight(length(rangeK),rangeK(end),c+(k-1)*rangeK(end),0.015)
            
            hold on
            
            LT_cond1 = squeeze(LT(Index_Conditions == condRow(i),rangeK == k+1,c))';  % Vector containing LT of c in cond1
            LT_cond2 = squeeze(LT(Index_Conditions == condCol(i),rangeK == k+1,c))';  % Vector containing LT of c in cond2
            
            if LT_pval2sided(i,k,c) > 0.05
                bar(1:2,[mean(LT_cond1) mean(LT_cond2)],'EdgeColor','k','LineWidth',0.8,'FaceColor','none')
            end
            if LT_pval2sided(i,k,c) < 0.05 && LT_pval2sided(i,k,c) > (0.05/rangeK(k))
                bar(1:2,[mean(LT_cond1) mean(LT_cond2)],'EdgeColor','r','LineWidth',0.8,'FaceColor','none')
            end
            if LT_pval2sided(i,k,c) < (0.05/rangeK(k)) && LT_pval2sided(i,k,c) > (0.05/sum(rangeK))
                bar(1:2,[mean(LT_cond1) mean(LT_cond2)],'EdgeColor','g','LineWidth',0.8,'FaceColor','none')
            end
            if LT_pval2sided(i,k,c) <= (0.05/sum(rangeK))
                bar(1:2,[mean(LT_cond1) mean(LT_cond2)],'EdgeColor','b','LineWidth',0.8,'FaceColor','none')
            end
        
            hold on
            errorbar([mean(LT_cond1) mean(LT_cond2)],[std(LT_cond1)/sqrt(numel(LT_cond1)) std(LT_cond2)/sqrt(numel(LT_cond2))],...,
                     'LineStyle','none','Color','k','CapSize',5,'LineWidth',0.4);
            if k == length(rangeK)
                set(gca,'XTick',1:2,'XTickLabel',{cond{condRow(i)} cond{condCol(i)}},'Fontsize',6,'TickLabelInterpreter','none')
                % xtickangle(45)
            else
                set(gca,'XTick',1:2,'XTickLabel',[],'Fontsize',6,'TickLabelInterpreter','none')
            end
            set(gca,'color','none')
        
            hold off
            box off
        end
    end
    saveas(Fig3, fullfile(data_dir, ['DwellTime_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.png']),'png');
    saveas(Fig3, fullfile(data_dir, ['DwellTime_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.fig']),'fig');
    disp(['- Plot successfully saved as DwellTime_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)}]);
end

% Plot of Hedge's effect size from dwell time hypothesis tests
disp(' ')
disp('Plotting Hedge''s effect size from hypothesis (permutation) tests:')
Fig4 = figure('Position', get(0, 'Screensize'));
for s_ind = 1:length(subplot_indices)
    
    subplot_ind = subplot_indices(s_ind);
    
    % Subplot adds plots following the rows (not the columns)
    subplot(size(subplot_map,1),size(subplot_map,2),index_fig(subplot_ind))
    
    plot(rangeK(1)-1:rangeK(end)+1,0.8*ones(1,length(rangeK)+2),'b--','LineWidth',1.5);
    hold on;
    plot(rangeK(1)-1:rangeK(end)+1,0.5*ones(1,length(rangeK)+2),'g--','LineWidth',1.5);
    plot(rangeK(1)-1:rangeK(end)+1,0.2*ones(1,length(rangeK)+2),'r--','LineWidth',1.5);
    
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            if effectsize(s_ind,k,c) < 0.2
                plot(rangeK(k),effectsize(s_ind,k,c),'.k','Markersize',4);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.2 && effectsize(s_ind,k,c) < 0.5
                plot(rangeK(k),effectsize(s_ind,k,c),'+r','Markersize',4);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.5 && effectsize(s_ind,k,c) < 0.8
                plot(rangeK(k),effectsize(s_ind,k,c),'og','Markersize',4);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
            if effectsize(s_ind,k,c) > 0.8
                plot(rangeK(k),effectsize(s_ind,k,c),'*b','Markersize',4);
                % text(rangeK(k),effectsize(s_ind,k,c),[' ' num2str(c)]);
            end
        end
    end
    
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}],'interpreter','none')
    ylabel('Hedge''s effect size')
    xlabel('Number of PL States K')
    xticks([2 4 6 8 10 12 14 16 18 20])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    ylim([0 1])
    set(gca,'YTick',[0,0.2,0.5,0.8,1])
end
saveas(Fig4, fullfile(data_dir, 'DwellTime_effetcsize.png'),'png');
saveas(Fig4, fullfile(data_dir, 'DwellTime_effetcsize.fig'),'fig');
disp('- Plot successfully saved as DwellTime_effetcsize');

disp(' ')
close all;
