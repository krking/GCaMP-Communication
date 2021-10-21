%% Crosscorrelation heatmap of rank-ordered calcium spikes for all cells

% binary metric based on active ROI and sort upon the most active ROI to
% least 
 m_cell_general = m_cell(cell_spike(2,:),:);
% signal metric based on active ROI and sort upon the most active ROI
 value_general = val_cell(cell_spike(2,:),:);
% Creating a new matrix to investigate  the first time of spike and the
% responsible roi
 [roi,tos,~] = find(m_cell_general ==1); 
 %Save both as a new array
 comb_t_roi = [ roi'; tos'] ;
 % Save the index of sorted rois from binary metric
 [~,ii] = sort(roi');
 comb_t_roi = comb_t_roi(:,ii);
 value_general = []; 
 for j = 1: length(comb_t_roi)
     value_general(:,j) = val_cell(cell_spike(2,comb_t_roi(1,j)), comb_t_roi(2,j)); 
 end
 %add the df/f value of each spiking to the new array
comb_t_roi = [comb_t_roi ; value_general]; 

%Spatial location of single cells with color indicating the timing of cellular calcium-dependent
%fluorescence spikes

c = jet(length(cell_spike));
roi = sort_cell(cell_spike(2,:));
figure(1);clf;scatter(s_loc(1,roi), s_loc(2,roi),114,c,'filled');
colormap(jet);
colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'1','4','7','10','13 spikes'})

%% Fig 3C 
%Color represents the product of spike time difference and Euclidean spatial distance

%sorting time of spiking 
[~,i] = sort(comb_t_roi(2,:));
time_relate= comb_t_roi(2,i); 
% importing all cells and sorting based on time of spiking 
cell_imp = sort_cell(cell_spike(2,comb_t_roi(1,:)));
roi_imp_real2 = cell_imp(i);
 s_loc_imp2 = s_loc(:, roi_imp_real2); 

 dist_imp_roi = pdist2(s_loc_imp2.',s_loc_imp2.');
 
 dist_time = pdist2(time_relate',time_relate');

 cluster_pair_diff =dist_time .* (dist_imp_roi.^1.5);

 figure(1);clf;
heatmap(cluster_pair_diff ,'Colormap',flipud(jet),'ColorScaling','log','ColorLimits',[0 15],'CellLabelColor','none');
 figure(2);clf;
heatmap(dist_time ,'Colormap',flipud(jet),'ColorScaling','log');
 figure(3);clf;
heatmap(dist_imp_roi ,'Colormap',flipud(jet));
