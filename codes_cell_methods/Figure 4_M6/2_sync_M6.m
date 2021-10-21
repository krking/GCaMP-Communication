%% Unsupervised k-mean clustering based on Euclidean distances  Figure 4B and C

% Finding ROIs that are sync together

% create S based on optimum window 
opt_w = 10; 
 t_bin = []; 
 
for  j = 1:(size(y2,2) - opt_w);
    [locs_pks, col] = find(m_cell(:,(j:j+opt_w)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    t_bin(:,j) = m ; 
          
end

[~,tseries] = find(t_bin>Sth(opt_w)); 
twindow = zeros(length(tseries), 2); 
twindow(:,1) = tseries  ; 
twindow(:,2) = tseries + opt_w; 


twindow(twindow(:,2)> 203,:) =[];
% which ROIs has spike within this range 
sync_roi = {};
for i = 1: length(twindow)
    t = twindow(i,:); 
    [roi, ~] = find(m_cell(:,t(1):t(2))==1); 
    sync_roi{i} = roi; 
end

sync_roi_vct = vertcat(sync_roi{:});
[sync_roi_vct,~,~] = unique(sync_roi_vct,'first');
sync_roi_vct2 = sort(sync_roi_vct);
sync_roi_vct = sort(sort_cell(sync_roi_vct)); 
sync_roi_str = string(sync_roi_vct);
sync_roi_loc= s_loc(:,sync_roi_vct.');

figure(1) ;clf; scatter(sync_roi_loc(1,:), sync_roi_loc(2,:),'green','filled');
dx = 0.3; dy = 0.3; % displacement so the text does not overlay the data points
text(sync_roi_loc(1,:)+dx,sync_roi_loc(2,:)+dy, sync_roi_str);
% xlim([ 0 636.16]);
% ylim([0 636.16]);
title 'centroid of synchronous cells' 

m_cell_wz_12 = m_cell(sync_roi_vct2,:); 

[~,b3] = find(comb_rnks_pks(2,:) ==sync_roi_vct);
% collect the number of spikes for sync ROIs
cb = comb_rnks_pks(:,b3);
%clustergram of sync ROI
cg2 = clustergram(y2(sync_roi_vct2,:),'RowLabels',sync_roi_str,'Colormap','jet',...
'Cluster',1,'Standardize',2);
% Figure 4B Unsupervised k-mean clustering
opts = statset('Display','final');
[idx,C,sumd,D] = kmeans(sync_roi_loc.',2,'Distance','sqeuclidean',...
    'Replicates',5);
a = sync_roi_loc.';


figure(2); clf; 
plot(a(idx==1,1),a(idx==1,2),'r.','MarkerSize',12)
hold on
plot(a(idx==2,1),a(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'k*',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
dx = 0.3; dy = 0.3; % displacement so the text does not overlay the data points
text(sync_roi_loc(1,:)+dx,sync_roi_loc(2,:)+dy, sync_roi_str);
hold off

%cluster1
cg3 = clustergram(y2(sync_roi_vct2(idx==1,:),:),'Colormap','jet',...
'Cluster',1,'Standardize',2);
%cluster2
cg4 = clustergram(y2(sync_roi_vct2(idx==2,:),:),'Colormap','jet',...
'Cluster',1,'Standardize',2);
