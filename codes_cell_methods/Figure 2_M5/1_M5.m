%% IDENTIFY PEAKS & SMOOTHING
%Reading CSV  file 
%change the pathway
% IMP IMP file :Results_15t_M07_set1
   
pathway = "/Applications/Calcium_file_Kings_lab/Results_013t_dish24_March2020.csv";
pathway2 = "/Applications/Calcium_file_Kings_lab/Results_bc_013t_dish24_march2020.csv"; 

bc_read = readtable(pathway2,'ReadVariableNames',true);
s_read = readtable(pathway,'ReadVariableNames',true);


% This table includes the centroid and mean signal
hasMatch = ~cellfun('isempty', regexp(s_read.Properties.VariableNames, 'Mean')) ;
s_read_short_mean = s_read(:,s_read.Properties.VariableNames(hasMatch));
s_read_location = s_read(:, s_read.Properties.VariableNames(~hasMatch));
s_read_short_location = s_read_location(1, 2:end);
xmatch =  ~cellfun('isempty', regexp(s_read_short_location.Properties.VariableNames, 'X')) ;
s_location_x = num2cell(  table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(xmatch))));
% check the movie and write down the size of the area x x micron
s_location_y = num2cell( 318.08 - table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(~xmatch))));
s = table2array(s_read_short_mean)';
s_loc = cell2mat(vertcat(s_location_x, s_location_y));

%background 
%bc for 013_dish24 
bc_read2 = table2array(bc_read(:,2))';

s = s - bc_read2; 
s(s<0) =0; 
s_pks_smooth2_imp = movmean(s,5,2); 
s_pks_smooth1 = movmean(s,200,2); 



%% plotting 
x= 1: size(s,2);

[signal_val,signal_to_pick] = max((s_pks_smooth2_imp - s_pks_smooth1),[],2);
conc_s = [signal_val,signal_to_pick];
[~,indx] = sort(signal_val);
conc_s = [conc_s(indx,:),indx] ;
topcell = conc_s(:,3); 
 c = topcell(40:50)';
 %c = [116 107 6 117];
cell_location = cell((size(c)));
figure(1); clf;
 for i = 1:length(c)
     cell_location(:,i) = {['Cell ', int2str(c(:,i))]};
    plot(x , s(c(:,i),:) + 100*(i-1)*ones(size(s(1,:))),...
    'LineWidth',1);
    hold on 
 end
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of CSF1r BMDM cells']);
hold off


%% IMP - filtering cells that are not spiking 
%remove 10% cells that are not active 
topcell = conc_s(:,3); 
topcell(1:round(size(s,1)/10)) = []; 
sort_cell = sort(topcell);

s_pks_smooth2_imp_2 =s(sort_cell,:);
figure(2);clf;
 h = heatmap(s_pks_smooth2_imp_2, 'ColorbarVisible', 'on', 'XLabel', 'Time', ...
    'YLabel', 'Cells', 'GridVisible','off','Colormap',jet, 'ColorScaling','scaledrows');
xlabel('Time Frame');
ylabel('Each individual ROIs');
title('DNA');
%now you can apply df/f and smooth it/ just remember if you wanna
%investigate the location, go back to sort_cell 
df_f = cal_df(s_pks_smooth2_imp_2);


%% denoising using zero filtering - low pass IIR
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.25,'DesignMethod','butter');
y2 = []; 

for i = 1: size(df_f,1)
y2(i,:) = filtfilt(d1,df_f(i,:));

end 

kk = 64 ;
k = find(sort_cell ==kk) 

figure(3);clf; subplot(3,1,1); plot(s(kk,:));  subplot(3,1,2); plot(df_f(k,:));...
subplot(3,1,3); plot(y2(k,:)); 

figure(4);clf; 
subplot(2,1,1)
plot(df_f(k,:))
hold on
plot(y2(k,:),'LineWidth',2)
legend('Noisy signal','Zero-Phase Filtering')
subplot(2,1,2) 
plot(s(kk,:)); 
legend('original signal'); 

%% finding peaks 
locs_zero_one = zeros(size(df_f));
width_zero_one = zeros(size(df_f));
prominence_zero_one = zeros(size(df_f));
m_cell2 = zeros(size(df_f)); 
pks_zero_val = zeros(size(df_f));

%varry parameter 
%CHANGING THEM 

th = 1e-4;

minpeakheight = 0.5;
for i = 1: size(df_f,1)
    for j = 1: size(df_f,2)
        [pks, locs, width, p] = findpeaks(y2(i,:),'Threshold', th,'MinPeakHeight', minpeakheight,'MinPeakDistance',10);
        m = size(locs,2);
        for n = 1:m
        if j == locs(:,n)
            locs_zero_one(i,j) = 1;
        end
        end
        for m2= 1:m
            if j ==locs(:,m2)
                width_zero_one(i,j) =width(:,m2); 
                 prominence_zero_one(i,j) =p(:,m2); 
                 pks_zero_val(i,j)= pks(:,m2);
                 m_cell2(i,j) =1; 
            end 
        end 
    end
  
end

prom_imp = prominence_zero_one; 

prom_imp(prom_imp ==0) = NaN; 

prom_imp = prom_imp(~isnan(prom_imp));
median(prom_imp)
mean(prom_imp)
m_cell = zeros(size(df_f));
val_cell= zeros(size(df_f));

width_new = zeros(size(df_f)); 
prominence_new = zeros(size(df_f)); 
for i = 1: size(df_f,2)
    for j = 1: size(df_f,1) 
        if (prominence_zero_one(j,i) > 0.5 && width_zero_one(j,i) <200)
            prominence_new(j,i) = prominence_zero_one(j,i); 
            width_new(j,i) = width_zero_one(j,i); 
            m_cell(j,i) =1; 
            val_cell(j,i) = y2(j,i);
        end
    end 
end 


num_pks = sum(m_cell.');
[rnk_pks,I] = sort(num_pks,'descend');
II = sort_cell(I); 
comb_rnks_pks = [rnk_pks ;II'];
comb_rnks_2 = [rnk_pks; I]; 
%which cell to check to investigate which ROIs are true positive 
k = I(2)
kk = II(2)


locs_new_peaks = find(m_cell(k,:) ==1)
figure(5); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk)]);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)

figure(6);clf; 
stem(m_cell(k,:),'g','filled');
ax = gca; 
ax.YLim = [0 1.5];
title('Inferred peaks');
figure(7); clf; 
histogram(rnk_pks);

[~,aa] = find(comb_rnks_pks(1,:) ~= 0);
wc_cm = sort(comb_rnks_pks(2,aa)); 
length(wc_cm)/ size(df_f,1)
sync_roi_str = string(wc_cm);
sync_roi_loc= s_loc(:,wc_cm.');
%display of active rois
figure(8) ;clf; scatter(sync_roi_loc(1,:), sync_roi_loc(2,:),'green','filled');
dx = 0.3; dy = 0.3; % displacement so the text does not overlay the data points
text(sync_roi_loc(1,:)+dx,sync_roi_loc(2,:)+dy, sync_roi_str);
% for 15t is 317.95
xlim([ 0 318.08]);
 ylim([0 318.08]);
title 'centroid of synchronous cells' 


%% histogram of time of spiking 
[row,col] = find(m_cell ==1);
a = [row, col]; 
[wcell,ind] = sort(row,'ascend');
tim_spk = col(ind);
comb_cor = [sort_cell(wcell); tim_spk];

 figure(9); clf; h2= histogram(tim_spk, 'BinWidth',50);
 
 edgesreal = 0:12.5:600; 
 nreal = histcounts(tim_spk,edgesreal);
 initiate_c = find(sort_cell ==41)
 k =initiate_c
kk = 41
locs_initiator = find(m_cell(k,:) ==1)
 s_initiate = y2(initiate_c,1:300); 
 %focus on cells that are spiking 
cell_spike = comb_rnks_2(:,find(rnk_pks~=0)); 
 %remove cell 12
sort_cell(10)
cell_spike(:,cell_spike(2,:) ==10) = [];

 m_cell_general = m_cell(cell_spike(2,:),:);
 value_general = val_cell(cell_spike(2,:),:);
 
 [roi,tos,~] = find(m_cell_general ==1); 
 
 [~,ii] = sort(roi');
 comb_t_roi = [ roi'; tos'] ;
 comb_t_roi = comb_t_roi(:,ii);
 value_general = []; 

 for j = 1: length(comb_t_roi)
     value_general(:,j) = val_cell(cell_spike(2,comb_t_roi(1,j)), comb_t_roi(2,j)); 
    
 end
comb_t_roi = [comb_t_roi ; value_general]; 

k = cell_spike(2,11)
kk = sort_cell(k)
locs_new_peaks = find(m_cell(k,:) ==1)
figure(10); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk)]);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)


maxdf = [];
maxdf4 =[]; 
ij= 0 ; 
for j= 1:length(cell_spike) 
    jj = find(comb_t_roi(1,:) ==j);
     maxdf4(:,j) = (0.95 *max(comb_t_roi(3,jj))); 
      [~, i] = min(abs(maxdf4(:,j)-comb_t_roi(3,jj)));
         maxdf(:,size(maxdf,2)+1) = comb_t_roi(:,i+ij);
         ij = jj(end);
end 
roi_spike = cell_spike(2, roi');
maxdfA = maxdf; 
maxdfA(1,:) = cell_spike(2, maxdf(1,:));
[~,maxdf3] = sort(maxdfA(2,:)); 

maxdfA= maxdfA(:,maxdf3); 
wc_ht = maxdfA(1,:); 
spk_ht = y2(wc_ht,:); 
%Figure 2C
figure(11);clf; 

h = heatmap(spk_ht, 'GridVisible','off','Colormap',jet,'ColorScaling','scaledrows');

 BIMP =maxdfA; 
  initiate_c = find(sort_cell ==41)

 find(BIMP(1,:) ==initiate_c)
 BIMP(2,find(BIMP(1,:) ==initiate_c)) = 118;
  %BIMP(:,BIMP(2,:) >150) = [];
 roi_imp_real = sort_cell(BIMP(1,:));
 
 s_loc_imp = s_loc(:,roi_imp_real);
 dist_roi = pdist2(s_loc_imp.',s_loc(:,41).'); 
 [~,i]= sort(BIMP(2,:));
 roi_imp_real2 = roi_imp_real(i);
time_relate= BIMP(2,i); 
 
[~,i] = sort(comb_t_roi(2,:));
time_relate= comb_t_roi(2,i); 
cell_imp = sort_cell(cell_spike(2,comb_t_roi(1,:)));
roi_imp_real2 = cell_imp(i);
 s_loc_imp2 = s_loc(:, roi_imp_real2); 
 %[318.08 318.08]
 %Figure2F
 dist_imp_roi = pdist2(s_loc_imp2.',s_loc_imp2.');
 
 dist_time = pdist2(time_relate',time_relate');

 cluster_pair_diff =dist_time .* (dist_imp_roi.^1.5);
[a,i] = find( cluster_pair_diff< 500 & cluster_pair_diff ~=0);
 figure(14);clf;
heatmap( cluster_pair_diff ,'Colormap',flipud(jet),'ColorScaling','log');
 figure(15);clf;
heatmap(dist_time ,'Colormap',flipud(jet),'ColorScaling','log');
 figure(16);clf;
heatmap(dist_imp_roi ,'Colormap',flipud(jet));