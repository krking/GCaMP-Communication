%% IDENTIFY PEAKS & SMOOTHING
%Reading CSV  file 
%change the pathway
% IMP IMP file :Results_15t_M07_set1  
pathway = "/Applications/Calcium_file_Kings_lab/Results_19t_M07_May2020.csv";
s_read = readtable(pathway,'ReadVariableNames',true);
pathway2 = "/Applications/Calcium_file_Kings_lab/Results_19t_M07_bc_detail_April2020.csv"; 
bc_read = readtable(pathway2,'ReadVariableNames',true);
% This table includes the centroid and mean signal
hasMatch = ~cellfun('isempty', regexp(s_read.Properties.VariableNames, 'Mean')) ;
s_read_short_mean = s_read(:,s_read.Properties.VariableNames(hasMatch));
s_read_location = s_read(:, s_read.Properties.VariableNames(~hasMatch));
s_read_short_location = s_read_location(1, 2:end);
xmatch =  ~cellfun('isempty', regexp(s_read_short_location.Properties.VariableNames, 'X')) ;
s_location_x = num2cell(  table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(xmatch))));
% check the movie and write down the size of the area x x micron
s_location_y = num2cell(317.95 -  table2array(s_read_short_location(:,s_read_short_location.Properties.VariableNames(~xmatch))));
s = table2array(s_read_short_mean)';
ds = median(s,2);
s_loc = cell2mat(vertcat(s_location_x, s_location_y));

%background 

 bc_read2 = table2array(bc_read(:,2))';

 bc_read2 = cal_df(bc_read2);
s = cal_df(s);
  df_f = s - bc_read2;
 df_f(df_f<0) = 0;

s_pks_smooth2_imp = movmean(df_f,5,2); 
s_pks_smooth1 = movmean(df_f,200,2); 

%% ploting 
%vector of time frame: 
[signal_val,signal_to_pick] = max((s_pks_smooth2_imp - s_pks_smooth1),[],2);
conc_s = [signal_val,signal_to_pick];
[~,indx] = sort(signal_val);
conc_s = [conc_s(indx,:),indx] ;
topcell = conc_s(:,3); 
x= 1:600;
%which cell to look 
c= topcell(518:523)';

cell_location = cell((size(c)));
figure(1); clf;
 for i = 1:length(c)
  cell_location(:,i) = {['Cell ', int2str(c(i))]};
    plot(x , df_f(c(i),:) + 0.5*(i-1)*ones(size(s(1,:))),...
    'LineWidth',1);
    hold on; 
 
 end
 
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of invivo cells']);
plot(bc_read2);
hold off;
%% removing ROIs that are not active 
% since all Active Rois have been collected, we attempted to remove only
% noisy signals
lowmax = conc_s(:,1)<= 0.2;
topcell = conc_s(lowmax ==0,3); 

sort_cell = sort(topcell);

df_f =df_f(sort_cell,:);
%Figure 5B 
figure(2);clf;
 h = heatmap(df_f, 'ColorbarVisible', 'on', 'XLabel', 'Time', ...
    'YLabel', 'Cells', 'GridVisible','off','Colormap',jet, 'ColorScaling','scaledrows');
xlabel('Time Frame');
ylabel('Each individual ROIs');
title('DNA');
%% denoising using zero filtering - low pass IIR 
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',1,'DesignMethod','butter');

y2 = []; 
for i = 1: size(df_f,1)
y2(i,:) = filtfilt(d1,df_f(i,:));

end 
y2(y2<0) =0 ;
kk = 41 ;
k = find(sort_cell ==kk) 

figure(3);clf; subplot(3,1,1); plot(s(kk,:));  subplot(3,1,2); plot(df_f(k,:));...
subplot(3,1,3); plot(y2(k,:)); 

figure(4);clf; 
subplot(2,1,1)
plot(df_f(k,:))
hold on
plot(y2(k,:),'LineWidth',1)
% hold on 
% plot(A2(k,:),'LineWidth',2,'Color','g');
legend('Noisy signal','Zero-Phase Filtering')
subplot(2,1,2) 
plot(s(kk,:)); 
legend('original signal'); 
%% findpeaks 
locs_zero_one = zeros(size(df_f));
width_zero_one = zeros(size(df_f));
prominence_zero_one = zeros(size(df_f));
pks_zero_val = zeros(size(df_f));
 m_cell2 = zeros(size(df_f)); 

%varry parameter 
%CHANGING THEM 
% for 15t the minpeakdist is 12 
th = 1e-4;
minpeakheight = 0.5;
minpeakdist = 5;
%for 15t_movie use 10 
%use s_pks_smooth2_imp
for i = 1: size(df_f,1)
    for j = 1: size(df_f,2)
        [pks, locs, width, p] = findpeaks(y2(i,:),'Threshold', th, 'MinPeakDistance',...
        minpeakdist,'MinPeakHeight', minpeakheight);
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
                % pks_zero_val(i,j)= pks(:,m2);
                 m_cell2(i,j) =1; 
            end 
        end 
    end
     
end

%findpeaks using parameters like size of width and prominence for 15t_
%movie not for 46t movie 

prom_imp = prominence_zero_one; 

prom_imp(prom_imp ==0) = NaN; 

prom_imp = prom_imp(~isnan(prom_imp));
mean(prom_imp)
median(prom_imp)
%or use 1 
par_w = 3.5; 
val_cell= zeros(size(df_f));
m_cell =zeros(size(df_f));

for i = 1: size(df_f,2)
    for j = 1: size(df_f,1)
if( width_zero_one(j,i) <50 && width_zero_one(j,i) >par_w && prominence_zero_one(j,i)>0.5)
 m_cell(j,i)  = 1;  
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

kk = find(II == 22)
kk =II(157)
k = I(157)



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

%% different method using number of peaks and time of peaks that happened 

%number of peak is: 
rnk_pks(rnk_pks ==0) =[];
figure(6); clf; 

histogram(rnk_pks, 19);

[row,col] = find(m_cell ==1);
a = [row, col]; 
[wcell,ind] = sort(row,'ascend');
tim_spk = col(ind);
comb_cor = [sort_cell(wcell); tim_spk];
 figure(7); clf; h2= histogram(tim_spk, 'BinWidth',30);

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

k = cell_spike(2,76)
kk = sort_cell(k)
locs_new_peaks = find(m_cell(k,:) ==1)
figure(8); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk)]);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)

A = [];
ij= 0 ; 
for j= 1:length(rnk_pks)
    jj = find(comb_t_roi(1,:) ==j);
      
     [~, i] = max(comb_t_roi(3,jj) - (mean(comb_t_roi(3,jj))- var(comb_t_roi(3,jj))));
        A(:,size(A,2)+1) = comb_t_roi(:,i+ij);
        ij = jj(end);
end 
roi_spike = cell_spike(2, roi');
AA = A; 
AA(1,:) = cell_spike(2, A(1,:));
[~,A3] = sort(AA(2,:)); 

AA= AA(:,A3); 
wc_ht = AA(1,:); 
spk_ht = m_cell(wc_ht,:); 

figure(9);clf;
h = heatmap(y2(wc_ht,:), 'GridVisible','off','Colormap',jet,'ColorLimits',[0 3]);

cg = clustergram(df_f,'AnnotPrecision',3,'Colormap','jet',...
'Cluster',1,'Standardize', 2);

pd_np = fitdist(rnk_pks.', 'NegativeBinomial');
[f, p, ss,cv] = kstest(rnk_pks.', 'CDF',pd_np)
othergrid_rnk = transpose(linspace(min(rnk_pks(:)),max(rnk_pks(:)),100));
figure(12);clf;
histogram(rnk_pks,13,'Normalization','pdf','FaceColor','g'); %plot original data
w_rnk = pdf(pd_np,othergrid_rnk);
hold on
plot(othergrid_rnk,w_rnk,'LineWidth',2,'Color','r') %plot GMM over original data
ylabel('Frequency'); 
xlabel('number of peaks'); 
title('Histogram of Number of peaks for each ROIs');
legend('Real data,number of peaks','NegativeBinomial fit');
hold off
% dont use truncate function for 15t movie 
pd_np= truncate(pd_np, 0,18);
rnd_np = random(pd_np,size(m_cell,1),1);
[nreal, edgesreal] = histcounts(rnk_pks,17);
figure(13); clf;
subplot(2,1,1); histogram(rnk_pks, 18,'FaceColor','g'); 
title('real number of peaks within each ROIs');
subplot(2,1,2); histogram(rnd_np,18,'FaceColor',[0.75 0.75 0.75]);
title('random number of peaks within each ROIs NegativeBinomial dist');
rnd_imp = zeros(size(m_cell,1),2000); 
n = 0 ; 
if n < size(m_cell,1)
for i =1 :2000
    j = 0;
j = random(pd_np,size(m_cell,1),1);
rnd_imp(:,i) = j; 

 
end
    [n, edges,bin] = histcounts(rnd_imp,17);
    n = round(n./2000);
end
if sum(n) < size(m_cell,1)
    a = size(m_cell_general,1) - sum(n);
    n(:,1)= n(:,1) +a; 
else 
    a = sum(n) - size(m_cell,1) 
    n(:,1) = n(:,1) -a ; 
end
bin_spks= [0:size(n,2)-1;n];
figure(14); clf;
subplot(2,1,1); histogram(rnk_pks, 'BinEdges',edges,'FaceColor','g'); 
title('real number of peaks within each ROIs');
subplot(2,1,2); h = histogram('BinEdges',edges,'BinCounts',n,'FaceColor',[0.75 0.75 0.75]);
title('random number of peaks within each ROIs NegativeBinomial dist');
% where the peak happened? 
[~,wpks] = find(m_cell==1);
figure(15);clf; h =histogram(wpks,'BinWidth',30);
% csvwrite('figure5_wpks.csv', wpks);
[a,bb,b] = histcounts(wpks,'BinMethod','auto');
figure(16);clf;
h4= histogram('BinCounts',a,'BinEdges',bb,'FaceColor',[0.75 0.75 0.75]);

%IMP 
%giving cells random number of spikes based on pdf of it 
 art_spk =[];
  h= 0; 
   for j = 1:size(bin_spks,2) 
      
       for m = h+1:h+bin_spks(2,j) 
          art_spk(m,:) = bin_spks(1,j);  
       end
       h = h + bin_spks(2,j);
   end

rng('shuffle');
rand_perm_spk = randperm(size(art_spk,1));
art_spk2 = art_spk; 
art_spk2(rand_perm_spk,:) = art_spk(:,:);

%new artificial binary matrix
k = 0;
art_cell = zeros(size(m_cell,1), length(m_cell),100);
for j = 1:100
for i = 1:(size(m_cell,1))
    k= 0 ;
    if (art_spk2(i,:) ~=0)
        k = art_spk2(i);
       rng('shuffle')
       rnd_numpks = randi([1 600],1,k);
        art_cell(i,rnd_numpks,j) = 1;
    end
end 
end 

[a,bb,b] = histcounts(wpks,'BinWidth',5);
[~,wpks_art] = find(art_cell==1);
[n, edges] = histcounts(wpks_art,300);
n = round(n./100); 
edges = round(edges./100);
figure(17);clf; histogram('BinCounts',n,'BinEdges',edges,'FaceColor', [0.75 0.75 0.75]); hold on ;
histogram('BinCounts',a,'BinEdges',bb,'FaceColor','g');

figure(18); clf;  histogram('BinCounts',a,'BinEdges',bb,'FaceColor','g'); hold on; ...
    plot(1:600, ones(1,600)*(mean(n)+std(a))); hold off;

%% synchronicity
%This will let you to find an optimal size of window 
per_art_cell = permute(art_cell,[2,1,3]);
num_pks_art = sum(per_art_cell);
C = reshape(num_pks_art,[],size(art_cell,1),1);
a = randi([1 100], 1,1);
num_pks_art2 = C(a,:);
[rnk_pks_art,I_art] = sort(num_pks_art2,'descend');

[~,wpks] = find(m_cell==1);
a = randi([1 100], 1,1)
[~,wpks_art] = find(art_cell ==1);
[n ,edges] = histcounts(wpks,12); 
[n2,edges2] = histcounts(wpks_art,12);
n2 = round(n2./100); 
edges2 = edges2./100; 

figure(29);clf;
histogram('BinCounts',n,'BinEdges',edges,'FaceColor','g');
hold on;
histogram('BinCounts',n2,'BinEdges',edges2,'FaceColor',[0.75 0.75 0.75]);
hold off;