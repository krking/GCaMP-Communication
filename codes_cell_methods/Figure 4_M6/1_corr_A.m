 %%  Correlation between location and timing of individual cells that spike during synchrony event
% construct a matrix that find a highest df/f for each active ROI and record the time of spike
 %row1: active roi, row2: time of spike in regarding to highest df/f, row3: highest df/f
maxdf = [];
ij= 0 ; 
for j= 1:comb_t_roi(1, end)  
    jj = find(comb_t_roi(1,:) ==j);
     [~, i] = max(comb_t_roi(3,jj) - (mean(comb_t_roi(3,jj))- var(comb_t_roi(3,jj))));
         maxdf(:,size(maxdf,2)+1) = comb_t_roi(:,i+ij);
         ij = jj(end);
end 

 % construct a matrix that find a highest df/f for each active ROI and
 % record the number of spikes for that ROI
 %row1:  roi, row2: number of spike , row3: highest df/f
maxdfA = maxdf; 
maxdfA(1,:) = cell_spike(2, maxdf(1,:));
%sort based on the df/f value 
[~,maxdf3] = sort(maxdfA(2,:)); 
maxdfA= maxdfA(:,maxdf3); 

  roi_imp_real  = sort_cell(maxdfA(1,:)); 
 roi_imp_real2  = sort_cell(maxdfA(1,:)); 
  s_loc_imp2 = s_loc(:, roi_imp_real2); 
dist_imp_roi = pdist2(s_loc_imp2.', [ 0 0]);
  time_relate = maxdfA(2,:) ;
  
  comb_difuse = [ dist_imp_roi, time_relate', roi_imp_real]; 
  
   figure(1) ;clf;
scatter(comb_difuse(:,2),comb_difuse(:,1),'green','filled'); hold on; 
dx = 0.3; dy = 0.3; 
 xlim([0 600]);
 