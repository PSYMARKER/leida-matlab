function [cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,centroids,rangeK)
%
% Function to calculate overlap between LEiDA K-means results obtained in
% Glasser, DBS80, AAL116 or AAL120 parcellation with the 7 Resting State
% Networks from Yeo et al. 2011
%
% INPUT:
% parcellation  parcellation used to segment the brain ('AAL116' or
%               'AAL120' or 'glasser378' or 'dbs80')
% n_areas       number of areas of the parcellation to consider
% centroids     centroids obtained from applying K-means
% rangeK        range of clustering solutions considered for analysis
%
% OUTPUT:
% cc_V_yeo7     correlation between Yeo nets and Centroids
% p_V_yeo7      p-value of each correlation coefficient
%
%  Note that the Yeo parcellation does not cover some subcortical areas 
%  nor the cerebellum (so LEiDA results obtained including these strucutres 
%  it may make sence to compute the overlap removing these areas.
%  i.e. in AALl116 compare only the first 90
%  i.e. in AAL120 compare only the first 94 
%  i.e. in glasser378 compare only the first 360 
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org


% 1) Define the Yeo Networks in the new Parcellation

% First load the mask of the chosen parcellation in MNI space 2mm
V_Parcels = struct2array(load(['ParcelsMNI2mm'],['V_' parcellation]));
V_Parcels((V_Parcels > n_areas)) = 0;

% Load the mask of the Yeo parcellation in MNI space 2mm
V_Yeo = struct2array(load('ParcelsMNI2mm','V_Yeo7'));


N_Yeo = max(V_Yeo(:)); % number of Yeo networks
Yeo_in_new_Parcel = zeros(N_Yeo,n_areas);

% Create 7 vectors representing the 7 Yeo RSNs in new Parcellation scheme
for n = 1:n_areas
    indn = V_Parcels == n;
    for Net = 1:7
        Yeo_in_new_Parcel(Net,n) = numel(find(V_Yeo(indn) == Net))/sum(indn(:));
    end
end
clear V_Yeo V_Parcels indn

% 2) Compare with the LEiDA results

% Vector to store the correlation coefficient
cc_V_yeo7=zeros(length(rangeK),max(rangeK),N_Yeo);
% Vector to store the p-value obtained from computing the correlation coef.
p_V_yeo7=ones(length(rangeK),max(rangeK),N_Yeo);

for k = 1:length(rangeK)

    % centroids = Kmeans_results (LEiDA_Clusters.m file)
    VLeida = centroids{k}.C(:,1:n_areas);
    
    for Centroid = 1:rangeK(k)
        
        % Here define what part of the Centroids to be compared
        V = VLeida(Centroid,:);
        
        % To set negative elements to zero
        V=V.*(V>0);
        
        for NetYeo = 1:N_Yeo
            % Consider only the positive values from cluster centroids
            % Correlation of 1st cluster centroid with any of the Yeo nets
            % is 0 (Global Mode)
            [cc, p]  = corrcoef(V(:),Yeo_in_new_Parcel(NetYeo,:));
  
            cc_V_yeo7(k,Centroid,NetYeo) = cc(2); % select 2nd entry of corr matrix
            p_V_yeo7(k,Centroid,NetYeo)  = p(2); % select 2nd entry of corr matrix
            
        end
    end
end


