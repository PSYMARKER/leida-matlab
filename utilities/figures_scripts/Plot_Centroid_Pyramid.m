function Plot_Centroid_Pyramid(data_dir,cond,parcellation,n_areas,cortex_dir)
%
% Plot of the rendered centroids on cortical space and overlap with Yeo
% functional brain networks (Yeo et al. 2011)
%
% INPUT:
% data_dir      directory where the centroids obtained for each clustering
%               solution are saved
% cond          tags of each condition considered in the experiment
% parcellation  parcellation used to segment the brain
% n_areas       number of areas of the parcellation to be considered
% cortex_dir    direction of observation of the rendering of the PL states
%
% OUTPUT:
% Fig           PL states rendered on a transparent cortex with their
%               corresponding overlap with the resting-state networks
%               defined by Yeo et al., (2011)
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with results for the dwell time (output from LEiDA_stats_DwellTime.m)
file_clusters = 'LEiDA_Clusters.mat';

% Load required data:
load([data_dir file_clusters],'Kmeans_results', 'rangeK');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% If n_Cond == 2 color the title of the glass brains according to FracOccup
file_FracOccup = 'LEiDA_Stats_FracOccup.mat';
if n_Cond == 2 % Load 2-sided pvals from FracOccup hypothesis tests
    if isfile([data_dir file_FracOccup])
        load([data_dir file_FracOccup], 'P_pval2sided');
        P_pval2sided = squeeze(P_pval2sided);
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT CENTROIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
disp('Computing overlap of centroids (PL states) with Yeo resting state networks')
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);

% Color of the Yeo networks as in the original paper (Yeo et al., 2011)
YeoColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation])); 
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Orientation of the cortex
disp(' ')
disp(['Rendering the centroids (using ' cortex_dir '):'])
Fig = figure('Position', get(0, 'Screensize'));
for k = 1:length(rangeK)
    
    % Matrix with the cluster centroids for a given K
    VLeida = Kmeans_results{k}.C;
    disp(['- K = ' num2str(rangeK(k))])
    
    for Centroid = 1:rangeK(k)
        
        % Vector of the coordinates of the cluster centroid
        V = VLeida(Centroid,:);
        V = round(V*1e2)/1e2;
        
        [cc_net,net] = max(cc_V_yeo7(k,Centroid,:));
        
        Centroid_Vol = zeros(size(Volume));
        
        for n = find(V >= 0)
            Centroid_Vol(Volume == n) = 1;
        end 
        sregion = smooth3(Centroid_Vol > 0);
        
        % View networks (use different sizes depending on number of
        % conditions)
        if isfile([data_dir file_FracOccup]) && n_Cond == 2
            subplot_tight(length(rangeK),rangeK(end),Centroid+(k-1)*rangeK(end),0.01)
        else
            subplot_tight(length(rangeK),rangeK(end),Centroid+(k-1)*rangeK(end),0.008)
        end
             
        hold on
        % First plot a transparent cortex
        cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
        reducepatch(cortexpatch,0.01);
        isonormals(scortex,cortexpatch);
        
        if p_V_yeo7(k,Centroid,net) < 0.05/rangeK(k) 
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
        else
            % Color in Black if no significant overlap with any of the Yeo resting-state networks
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %'FaceAlpha', V(n))
        end
        
        material dull
        lighting gouraud
        switch cortex_dir
            case 'TopView'
                view(-90,90)
            case 'SideView'
                view(0,0)
        end
        
        daspect([1 1 1])
        camlight;
        xlim([5 105])
        ylim([5 85])
        zlim([10 80])
        axis off
        
        if isfile([data_dir file_FracOccup]) && n_Cond == 2 % if only 2 conditions color the title of each brain according to pval from FracOccup hyp. tests
            if P_pval2sided(k,Centroid) > 0.05
                title(['p<' num2str(P_pval2sided(k,Centroid)*10,'%1.0e')],'Fontsize',5,'Color','k')
            elseif P_pval2sided(k,Centroid) < 0.05 && P_pval2sided(k,Centroid) > (0.05/rangeK(k))
                title(['p<' num2str(P_pval2sided(k,Centroid)*10,'%1.0e')],'Fontsize',5,'Color','r')
            elseif P_pval2sided(k,Centroid) < 0.05/rangeK(k) && P_pval2sided(k,Centroid) > (0.05/sum(rangeK))
                title(['p<' num2str(P_pval2sided(k,Centroid)*10,'%1.0e')],'Fontsize',5,'Color','g')
            elseif P_pval2sided(k,Centroid) <= (0.05/sum(rangeK))
                title(['p<' num2str(P_pval2sided(k,Centroid)*10,'%1.0e')],'Fontsize',5,'Color','b')
            end
        end
    end
end
saveas(Fig, fullfile(data_dir, ['Centroid_Pyramid_' cortex_dir '.png']),'png');
saveas(Fig, fullfile(data_dir, ['Centroid_Pyramid_' cortex_dir '.fig']),'fig');
disp(['- Plot successfully saved as Centroid_Pyramid_' cortex_dir]);

disp(' ')
close all;
        
