%%IDENTIFY PEAKS & SMOOTHING
%Reading CSV file 
%using same ROIS that adapts to DNA epi data
pathway = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postDNA10-131_Stage1.csv";
%background of DNA stimulated 
pathway2 = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postDNA10-131_Stage1_bc.csv";
pathway3 = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postgcamp1_Stage1.csv";
%background of nonDNA stimulated 
pathway4 = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postgcamp1_Stage1_Area.csv";
pathway5 = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postDNA10-131_Stage1_Area.csv";
pathway6 = "/Applications/Calcium_file_Kings_lab/Epi/Epi_Csf1r_Gcamp5_postgcamp1_Stage1_bc.csv";
%diff_DIC
pathway7 = "/Applications/Calcium_file_Kings_lab/Epi/Results_diff_DIC_Stage1.csv";
dic_read =readtable(pathway7,  'ReadVariableNames',true);

%pathway = "/Users/kingslab/M07_invivo/Results_19t_M07.csv"
s_read = readtable(pathway,'ReadVariableNames',true);
bc_read_dna =readtable(pathway2,'ReadVariableNames',true);
bc_read_dna2 = table2array(bc_read_dna(:,3));
s_read3 = readtable(pathway3,'ReadVariableNames',true);
bc_read_nodna = readtable(pathway6, 'ReadVariableNames',true);
bc_read_nodna2 = table2array(bc_read_nodna(:,3));
Area_whole = table2array(bc_read_nodna(1,2));
a_read = readtable(pathway5,'ReadVariableNames',true);
a2_read = readtable(pathway4,'ReadVariableNames',true);
Area_dna = table2array(a_read(:,2));
Area_nodna = table2array(a2_read(:,2));
%a(128,:) is the whole field 
%% DNA data
% This table includes the centroid and mean signal
hasMatch = ~cellfun('isempty', regexp(s_read.Properties.VariableNames, 'Mean')) ;
s_read_short_mean = s_read(:,s_read.Properties.VariableNames(hasMatch));
s_read_location = s_read(:, s_read.Properties.VariableNames(~hasMatch));
s_read_short_location = s_read_location(1, 2:end);
xmatch =  ~cellfun('isempty', regexp(s_read_short_location.Properties.VariableNames, 'X')) ;
s_location_x = table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(xmatch)));
s_location_y = table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(~xmatch)));
s = table2array(s_read_short_mean)';
s_loc = [s_location_x; s_location_y];
%distance in respect to ROI 3
distance = sqrt( (s_loc(1,:)- s_loc(1,3)).^2 + (s_loc(2,:)-s_loc(2,3)).^2 );
x = 1:size(s,1);
distance = [x; distance];
s = cal_df(s);
bc_read_dna2 = cal_df(bc_read_dna2');

cs_pks =s - bc_read_dna2;
cs_pks(cs_pks <0 )= 0;
cs_pks_norm = cs_pks./max(cs_pks,[],2);
s_pks_smooth2_imp = movmean(cs_pks,5,2); 
s_pks_smooth1 = movmean(cs_pks,200,2); 
%why removing 128 and 129? 
distance(:,128:129) = [];

%cs_pks= cs_pks(1:127,:);
hasMatch_dic = ~cellfun('isempty', regexp(dic_read.Properties.VariableNames, 'Mean')) ;
diff_dic = dic_read(:,dic_read.Properties.VariableNames(hasMatch_dic));
diff_dic = table2array(diff_dic)';
diff_dic_norm = diff_dic./ max(diff_dic,[],2);

dic_smooth2 = movmean(diff_dic,5,2); 
dic_smooth1 = movmean(diff_dic,501,2);
dic_imp = diff_dic - dic_smooth1;
%%
x= 1: size(cs_pks,2);
 c = [50];


%c = [116 107 6 117];
cell_location = cell((size(c)));
figure(1); clf;
 for i = 1:length(c)
     cell_location(:,i) = {['Cell ', int2str(c(:,i))]};
    plot(x , cs_pks(c(:,i),:) + 0.25*(i-1)*ones(size(cs_pks(1,:))),...
    'LineWidth',2);
    hold on 
 end
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of CSF1r BMDM cells']);
hold off

%% IMP - filtering cells that are not spiking 


[signal_val,signal_to_pick] = max((s_pks_smooth2_imp - s_pks_smooth1),[],2);
conc_s = [signal_val,signal_to_pick];
[~,indx] = sort(signal_val);
conc_s = [conc_s(indx,:),indx] ;

lowmax = conc_s(:,1)<= 0.02;

topcell = conc_s(:,3); 
topcell(find(lowmax ==1)) = []; 
sort_cell = sort(topcell);

s_pks_smooth2_imp_2 =cs_pks(sort_cell,:);
%% find peaks 
  %varry parameter 
locs_zero_one = zeros(size(s_pks_smooth2_imp_2));
width_zero_one = zeros(size(s_pks_smooth2_imp_2));
prominence_zero_one = zeros(size(s_pks_smooth2_imp_2));
features = cell(size(s_pks_smooth2_imp_2,1),1);
pks_zero_val = zeros(size(s_pks_smooth2_imp_2));
m_cell2 = zeros(size(s_pks_smooth2_imp_2));
%CHANGING THEM 
th = 1e-4; 
%use s_pks_smooth2_imp
for i = 1: size(s_pks_smooth2_imp_2,1)
    for j = 1: size(s_pks_smooth2_imp_2,2)
        [pks, locs, width, p] = findpeaks(s_pks_smooth2_imp_2(i,:),'Threshold', th, 'MinPeakDistance',5,'MinPeakHeight',0.01);
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
     features{i} = [width_zero_one(i,:); prominence_zero_one(i,:);pks_zero_val(i,:)];
end

%findpeaks using parameters like size of width and prominence 
par_w_min = 1; 
par_w_max = 10; 
prom_imp = prominence_zero_one; 

prom_imp(prom_imp ==0) = NaN; 

prom_imp = prom_imp(~isnan(prom_imp));
mean(prom_imp)

par_prom= 8;
m_cell = zeros( size(s_pks_smooth2_imp_2));
val_cell= zeros(size(s_pks_smooth2_imp_2));

for i = 1: size(s_pks_smooth2_imp_2,2)
    for j = 1: size(s_pks_smooth2_imp_2,1)
if( width_zero_one(j,i) >= par_w_min && prominence_zero_one(j,i)>0.049)
 m_cell(j,i)  = 1;  
 val_cell(j,i) = s_pks_smooth2_imp(j,i);
 width_new(j,i) = width_zero_one(j,i);
end 
    end
end
val_cell_norm = val_cell./ max(val_cell,[],2);
%% FIGURE 1 C &D
% total number of peaks for each cell

num_pks = sum(m_cell.');
[rnk_pks,I] = sort(num_pks,'descend');
II = sort_cell(I); 
comb_rnks_pks = [rnk_pks ;II'];
comb_rnks_2 = [rnk_pks; I]; 

 sz_bin = 1; 
 t_bin = []; 
 
for  j = 1:(size(s_pks_smooth2_imp_2,2) - sz_bin);
    [locs_pks, col] = find(m_cell(:,(j:j+sz_bin)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    t_bin(:,j) = m ; 
          
end
 t_bin = t_bin'; 

 figure(2);clf; bar(t_bin);
% to calculate the gfp area we need to figure it out which cells to
% consider
[~,col_c] = find(rnk_pks~=0);

gfp_k = I;
% which sorted cell
gfp_k = gfp_k(:,col_c);
% find which cell 
gfp_c= sort_cell(gfp_k);

% area gfp + FIGURE 1C

Area_DNA_gfp = (sum(Area_dna(gfp_c,:))/Area_whole)*100

% how many cells based on percentage FIGURE 1 D
(size(gfp_c,1)/ size(cs_pks,1))*100
 % to visualize peaks 
 
k = comb_rnks_2(2,10)
kk = sort_cell(k)
locs_new_peaks = find(m_cell(k,:) ==1)

figure(3); clf; 
plot(s_pks_smooth2_imp_2(k,:));hold on;
plot(locs_new_peaks, s_pks_smooth2_imp_2(k,locs_new_peaks),'r*'); hold on;

ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk),' of invivo studies']);
legend('calcium signal', 'peaks'); hold off;

prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)

%% max FIGURE 1E &H
cell_spike = comb_rnks_2(:,find(rnk_pks~=0)); 

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
maxdf = [];
maxdf4 =[]; 
ij= 0 ; 
for j= 1:length(cell_spike) 
    jj = find(comb_t_roi(1,:) ==j);
     [~, i] = max(comb_t_roi(3,jj) - (mean(comb_t_roi(3,jj)) -var(comb_t_roi(3,jj))));
         maxdf(:,size(maxdf,2)+1) = comb_t_roi(:,i+ij);
         ij = jj(end);
end 
roi_spike = cell_spike(2, roi');
maxdfA = maxdf; 
maxdfA(1,:) = cell_spike(2, maxdf(1,:));
[~,maxdf3] = sort(maxdfA(2,:)); 
maxdfA= maxdfA(:,maxdf3); 
wc_ht = maxdfA(1,:); 
spk_ht = s_pks_smooth2_imp_2(wc_ht,:); 

%FIGURE 1-E
figure(4);clf; 
h = heatmap(spk_ht, 'GridVisible','off','Colormap',jet,'ColorScaling','scaledrows','ColorLimits',[-0.3 1]);

%FIGURE 1-H
c = jet(length(wc_ht));
a  = sort_cell(wc_ht);
figure(5);clf;scatter(s_loc(1,sort_cell(wc_ht)), s_loc(2,sort_cell(wc_ht)),200,c,'filled','MarkerEdgeColor','black','LineWidth',1);
colormap(jet);
colorbar('Ticks',[0,0.1,0.2,0.3,0.4, 0.5, 0.6,0.7,0.8,0.9, 1],...
         'TickLabels',{'0','100','200','300','400','500','600','700','800','900','1000 min'})
%% difference between diff DIC and GFP FIGURE 1J 
[time_r, cell_r] = find(m_cell' ==1) ;
cell_r = sort_cell(cell_r);
[~,cell_r2,~] = unique(cell_r);
%first strategy 
time_1 = time_r - 19; 
time_2 = time_r +19; 
comb_import = [cell_r time_r time_1 time_2]; 
comb_import = comb_import(cell_r2,:);
%second strategy 
[~,i] = sort(maxdfA(1,:));
cell_imp = maxdfA(1,i)';
time_imp = maxdfA(2,i)';
time_1 = time_imp - 19; 
time_2 = time_imp +19; 
comb_import = [cell_imp time_imp time_1 time_2]; 


diff_dic_2 = zeros(length(comb_import),39); 
cs_pks_2 = zeros(length(comb_import),39);
for i = 1:length(comb_import)
 if comb_import(i,4) < 500 && comb_import(i,3) > 0 
diff_dic_2(i,:) = diff_dic(comb_import(i,1),comb_import(i,3):comb_import(i,4));

cs_pks_2(i,:) = cs_pks(comb_import(i,1),comb_import(i,3):comb_import(i,4));
 end
 
end
diff_dic_3 = diff_dic_2 ./ max(diff_dic_2,[],2);
cs_pks_norm = cs_pks_2./max(cs_pks_2,[],2);
wc = 10; 
figure(6);clf;subplot(2,1,1); plot(diff_dic_3(wc,:)); hold on;subplot(2,1,2); plot(cs_pks_norm(wc,:));
[~,pre_death] = max(diff_dic_2, [],2);
[~,post_death] = min(diff_dic_2,[],2);
[~,time_spk] = max(cs_pks_2,[],2);


b2 = diff_dic_3;
b3 = isnan(b2);
b3 = b3(:,1);
b2(b3 ==1,:) =[];
y2 = b2 -  min(mean(b2,1));
y3 = max(mean(b2,1)) - min(mean(b2,1));
y = y2./ y3;
csvwrite('figure1_1_mean_std_dic_6.csv',y);
csvwrite('figure1_1_mean_std_gfp_6.csv',cs_pks_norm);


%% noDNA data 

hasMatch3 = ~cellfun('isempty', regexp(s_read3.Properties.VariableNames, 'Mean')) ;
s_read_short_mean3 = s_read3(:,s_read3.Properties.VariableNames(hasMatch3));
s_read_location_3 = s_read3(:, s_read3.Properties.VariableNames(~hasMatch3));
s_read_short_location_3 = s_read_location_3(1, 2:end);
xmatch3 =  ~cellfun('isempty', regexp(s_read_location_3.Properties.VariableNames, 'X')) ;
s_location_x_3 = s_read_location_3(:,s_read_location_3.Properties.VariableNames(xmatch3));
s_short_location_x_3 = table2array(s_location_x_3(1,:));
s_location_y_3 = s_read_location_3(:,s_read_location_3.Properties.VariableNames(~xmatch3));
s_short_location_y_3 = table2array(s_location_y_3(1,2:end));
s3 = table2array(s_read_short_mean3)';
s_loc_3 = table2array(s_read_short_location_3);
s_min_3 = min(s3,[],2);

s3 = cal_df(s3);
bc_read_nodna2 = cal_df(bc_read_nodna2');

s_pks_3 = s3 - bc_read_nodna2;
s_pks_3(s_pks_3 <0 )= 0;
s_pks_smooth2_imp_nd = movmean(s_pks_3,5,2); 
s_pks_smooth1_nd = movmean(s_pks_3,200,2); 
figure(7);clf; plot(s3(46,:)); hold on ; plot(bc_read_nodna2'); hold on;...
plot(s_pks_3(46,:)); legend('s3','bc', 'smooth');

%% plotting nd data 
x= 1: size(s_pks_3,2);
 c = [51 46 5 63];
cell_location = cell((size(c)));
figure(8); clf;
 for i = 1:length(c)
     cell_location(:,i) = {['Cell ', int2str(c(:,i))]};
    plot(x , s_pks_3(c(:,i),:) + 0.2*(i-1)*ones(size(s_pks_3(1,:))),...
    'LineWidth',2);
    hold on 
 end
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of CSF1r BMDM cells']);
hold off

%% IMP - filtering cells that are not spiking 

[signal_valnd,signal_to_pick_nd] = max((s_pks_smooth2_imp_nd - s_pks_smooth1_nd),[],2);
conc_s2 = [signal_valnd,signal_to_pick_nd];
[~,indx_nd] = sort(signal_valnd);
conc_s2 = [conc_s2(indx_nd,:),indx_nd] ;

lowmax2 = conc_s2(:,1)<= 0.02;

topcell2 = conc_s2(:,3); 
topcell2(find(lowmax2 ==1)) = []; 
sort_cell2 = sort(topcell2);

s_pks_smooth2_imp_3 =s_pks_3(sort_cell2,:);


%% finding peaks for nd 
locs_zero_one_nd = zeros(size(s_pks_smooth2_imp_3));
width_zero_one_nd = zeros(size(s_pks_smooth2_imp_3));
prominence_zero_one_nd = zeros(size(s_pks_smooth2_imp_3));
features_nd = cell(size(s_pks_smooth2_imp_3,1),1);
pks_zero_val_nd = zeros(size(s_pks_smooth2_imp_3));
%CHANGING THEM 
th = 1e-4; 
%use s_pks_smooth2_imp
for i = 1: size(s_pks_smooth2_imp_3,1)
    for j = 1: size(s_pks_smooth2_imp_3,2)
        [pks, locs, width, p] = findpeaks(s_pks_smooth2_imp_3(i,:),'Threshold', th, 'MinPeakDistance',5,'MinPeakHeight',0.01);
        m = size(locs,2);
        for n = 1:m
        if j == locs(:,n)
            locs_zero_one_nd(i,j) = 1;
        end
        end
        for m2= 1:m
            if j ==locs(:,m2)
                width_zero_one_nd(i,j) =width(:,m2); 
                 prominence_zero_one_nd(i,j) =p(:,m2); 
                 pks_zero_val_nd(i,j)= pks(:,m2);
                 
            end 
        end 
    end
     features_nd{i} = [width_zero_one_nd(i,:); prominence_zero_one_nd(i,:);pks_zero_val_nd(i,:)];
end

%findpeaks using parameters like size of width and prominence 

m_cell3 = zeros( size(s_pks_smooth2_imp_3));
val_cell3= zeros(size(s_pks_smooth2_imp_3));
for i = 1: size(s_pks_smooth2_imp_3,2)
    for j = 1: size(s_pks_smooth2_imp_3,1)
if( width_zero_one_nd(j,i) >par_w_min  && prominence_zero_one_nd(j,i)>0.049)
 m_cell3(j,i)  = 1;  
 val_cell3(j,i) = s_pks_smooth2_imp_3(j,i);
end 
    end
end

num_pks3 = sum(m_cell3.');
[rnk_pks3,I3] = sort(num_pks3,'descend');
comb_rnks_pks3 = [rnk_pks3 ;I3];

[~,col_c_nd] = find(rnk_pks3~=0);

gfp_k_nd = I3;
gfp_k_nd = gfp_k_nd(:,col_c_nd);
gfp_c_nd= sort_cell2(gfp_k_nd);

%FIGURE 1C&D
Area_NODNA_gfp = (sum(Area_nodna(gfp_c_nd,:))/Area_whole)*100
Cell_NODNA_gfp = (length(gfp_c_nd)/size(s_pks_3,1))*100
k = comb_rnks_pks3(2,2)

kk = sort_cell2(k)
locs_new_peaks_nd = find(m_cell3(k,:) ==1);

figure(9); clf; 
plot(s_pks_smooth2_imp_3(k,:));hold on;
plot(locs_new_peaks_nd, s_pks_smooth2_imp_3(k,locs_new_peaks_nd),'r*'); hold on;

ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk),' of invivo studies']);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one_nd(k,locs_new_peaks_nd)
width_zero_one_nd (k,locs_new_peaks_nd)


%% which cell to plot for epi data/paper 1G
x = 1: length(s);
[a,b] = find(m_cell(106,:) ==1)

c= [3,4,6,14,46,47,63,65,17,88,95,97,109];



cell_location = cell((size(c)));
s_imp = zeros(length(c), length(cs_pks));

 for i = 1:length(c)
    s_imp(i,:) = s_pks_smooth2_imp_2(c(:,i),:) + 0.12*(i-1)*ones(size(cs_pks(1,:)));

 end
 figure(10);clf; plot(x, s_imp); 

csvwrite('firgure1_spikes.csv', s_imp');    