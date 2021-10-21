%% findpeaks 
locs_zero_one = zeros(size(df_f));
width_zero_one = zeros(size(df_f));
prominence_zero_one = zeros(size(df_f));
features = cell(size(df_f,1),1);
pks_zero_val = zeros(size(df_f));
m_cell2 = zeros(size(df_f)); 

%varry parameter_ to understand bette variation, please check findpeaks
%function from matlab 
%CHANGING THEM 
% for 15t the minpeakdist is 12 
th = 1e-5;
minpeakheight = 0.5;
minpeakdist = 5;
%This will take about 10-15 min to process 
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
                 pks_zero_val(i,j)= pks(:,m2);
                 m_cell2(i,j) =1; 
            end 
        end 
    end
  
end

% based on width and prominence build a binary signal matrix
val_cell= zeros(size(df_f));
m_cell =zeros(size(df_f));
par_w = 3.5;
for i = 1: size(df_f,2)
    for j = 1: size(df_f,1)
if( width_zero_one(j,i) <50 && width_zero_one(j,i) >par_w && prominence_zero_one(j,i)>0.5)
 m_cell(j,i)  = 1;  
 val_cell(j,i) = y2(j,i);
end 
    end
end

% to investigate which ROIs spike and how many times? 
num_pks = sum(m_cell.');
[rnk_pks,I] = sort(num_pks,'descend');
II = sort_cell(I); 
% Final array based on ROI number 
comb_rnks_pks = [rnk_pks ;II'];
% array based on ROIs that was above background
comb_rnks_2 = [rnk_pks; I]; 

%which cell to check to investigate which ROIs are true positive. ROI 197
k2 = find(II ==197)

kk =II(k2)
k = I(k2)
%Quality control based on visual plotting and checking movie
% at what time(s) spike
locs_new_peaks = find(m_cell(k,:) ==1)
figure(1); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk)]);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)

%% different method using number of peaks and time of peaks that happened 

%number of peak is:
zero_cell = size(find(rnk_pks ==0));
figure(2); clf; 
histogram(rnk_pks)


[row,col] = find(m_cell ==1);
[wcell,ind] = sort(row,'ascend');
%time of spiking from all active ROIs
tim_spk = col(ind);
% time of spiking + number of ROI that is involved
comb_cor = [sort_cell(wcell) tim_spk];
figure(3); clf; h2= histogram(tim_spk, 'BinWidth',30);
% number of spikes + number of active ROI
cell_spike = comb_rnks_2(:,find(rnk_pks~=0)); 
%number of active ROI
k = cell_spike(2,68)
% number of real ROI
kk = sort_cell(k)
locs_new_peaks = find(m_cell(k,:) ==1)
figure(4); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
% set(gca,'ytick',[])
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(kk)]);
legend('calcium signal', 'peaks'); hold off;
prominence_zero_one(k,locs_new_peaks)
width_zero_one(k,locs_new_peaks)


