%% IDENTIFY PEAKS & SMOOTHING
%Reading CSV  file 
%change the pathway

%mean value of ROI based on time
pathway = "/Applications/Calcium_file_Kings_lab/pathways/M6_meanvalue.csv";
s_read = readtable(pathway,'ReadVariableNames',true);
sread_short = s_read(:,2:end); 
%convert to array 
s = table2array(sread_short);
%remove noise sutble background 
pathway2 = "/Applications/Calcium_file_Kings_lab/pathways/M6_Background.csv"; 
bc_read = readtable(pathway2,'ReadVariableNames',true);
%convert to array 
bc_read2 = table2array(bc_read(:,2))';
%location 
p_s = "/Applications/Calcium_file_Kings_lab/pathways/M6_centroid.csv"  
p_sread = readtable(p_s, 'ReadVariableNames',true);
p_sread_short = p_sread(1,2:end); 
xmatch =  ~cellfun('isempty', regexp(p_sread_short.Properties.VariableNames, 'X')) ;

%location based on x and y and reconstruct the array
s_location_x = num2cell(  table2array(p_sread_short(:,p_sread_short.Properties.VariableNames(xmatch))));

s_location_y = num2cell(254.44 -  table2array(p_sread_short(:,p_sread_short.Properties.VariableNames(~xmatch))));

s_loc = cell2mat(vertcat(s_location_x, s_location_y));

figure(1);clf; subplot(1,2,1);plot(s(197,:)); hold on; ...
   plot(bc_read2);subplot(1,2,2);plot(s(321,:)); hold on; ...
   plot(bc_read2)

%background and signal in terms of df/f , remove the sutble background here
bc_read2 = cal_df(bc_read2);
s = cal_df(s);
df_f = s - bc_read2;
% remove the noise
df_f(df_f<0) = 0;

% smoothing signal using moving average filter 
%fine filter
s_pks_smooth2_imp = movmean(df_f,5,2); 
%coarse filter
s_pks_smooth1 = movmean(df_f,200,2); 
%% removing ROIs that are not active 
[signal_val,signal_to_pick] = max((s_pks_smooth2_imp - s_pks_smooth1),[],2);
conc_s = [signal_val,signal_to_pick];
[~,indx] = sort(signal_val);
conc_s = [conc_s(indx,:),indx] ;

topcell = indx; 
%remove top10%
prc = round(size(s,1)/10);
topcell(1:prc) = []; 

sort_cell = sort(topcell);

% ploting 
%vector of time frame: 
t= 1:600;

c = topcell(100:106,:)';
cell_location = cell((size(c)));
figure(1); clf;
 for i = 1:length(c)
     cell_location(:,i) = {['Cell ', int2str(c(:,i))]};
    plot(t , df_f(c(:,i),1:600) + 1.5*(i-1)*ones(size(s(1,1:600))),...
    'LineWidth',1.5);
    hold on; 
 
 end
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of invivo cells']);
hold off;

%removing ROIs 
df_f =df_f(sort_cell,:);

%heatmap 
figure(2);clf;
 h = heatmap(df_f, 'ColorbarVisible', 'on', 'XLabel', 'Time', ...
    'YLabel', 'Cells', 'GridVisible','off','Colormap',jet, 'ColorScaling','scaledrows');
xlabel('Time Frame');
ylabel('Each individual ROIs');
title('DNA');
%% denoising using zero filtering low pass IIR
d1 =designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.45,'DesignMethod','butter');
y2 = []; 

for i = 1: size(df_f,1)
y2(i,:) = filtfilt(d1,df_f(i,:));

end 
%plotting| ROI 187 in M6 video. Finding the respected roi in denoised
%signal 
kk = 187 ;
k = find(sort_cell ==kk) 

figure(3);clf; subplot(3,1,1); plot(s(kk,:));  subplot(3,1,2); plot(df_f(k,:));...
subplot(3,1,3); plot(y2(k,:)); 

figure(4);clf; 
subplot(2,1,1)
plot(df_f(k,:))
hold on
plot(y2(k,:),'LineWidth',3)
% hold on 
% plot(A2(k,:),'LineWidth',2,'Color','g');
legend('Noisy signal','Zero-Phase Filtering')
subplot(2,1,2) 
plot(s(kk,:)); 
legend('original signal'); 

%% clustergram 

cg = clustergram(df_f,'AnnotPrecision',3,'Colormap','jet',...
'Cluster',1,'Standardize', 2);
% to investigate more about the order of ROIs
sync_roi = cg.RowLabels;
sync_roi = sort_cell(str2double(sync_roi));
