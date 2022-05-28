%% clear your workspace
clc; 
clear all; 
%% Reading CSV  file 
%mean value of  each ROI raw fluorescence values for each frame 

[file,path]=uigetfile('*.csv') 
pathway = append(path,file)
%costumize the pathway based on the location of your storage


% read mean raw fluorescence values
s_read = readtable(pathway,'ReadVariableNames',true);
 Index_m = find(contains(s_read.Properties.VariableNames,'Mean'));
sread_short = s_read(:,Index_m); 
%convert to array 
s = table2array(sread_short).';
% read X and Y coordinate 
Index_x = find(contains(s_read.Properties.VariableNames,'X'));
s_location_x = table2array(s_read(1,Index_x));
Index_y = find(contains(s_read.Properties.VariableNames,'Y'));
s_location_y = table2array(s_read(1,Index_y));
% reconstructed X and Y array
s_loc = vertcat(s_location_x, s_location_y);
%%  Step 1: 
% %%constructing df/f: 
% bc_read2 = cal_df(bc_read2);
med_s =  median(s,2);                 
df_f = (s-med_s)./med_s;
% df_f = s - bc_read2;
% remove the background
df_f(df_f<0) = 0;

%%removing ROIs that are not active 

%fine filter
s_pks_smooth2_imp = movmean(df_f,5,2); 
%coarse filter
s_pks_smooth1 = movmean(df_f,200,2); 
 % ordering cells based on activity ( low to high) 
[signal_val,signal_to_pick] = max((s_pks_smooth2_imp - s_pks_smooth1),[],2);
conc_s = [signal_val,signal_to_pick];
[~,indx] = sort(signal_val);
conc_s = [conc_s(indx,:),indx] ;

topcell = indx; 
% sorting cells 
sort_cell = sort(topcell);
% ploting as an example 
%vector of time frame: 
t= 1:size(df_f,2);
%%INPUT for figure, which cells to display, the higher number the more
%%active cells are LINE 66
disp('We ordered ROIs based on activity and we represent based on index. The higher you go the more active ROIs are, Example: 110:118')
prompt = {'active index number of ROI 1:','active index number of ROI 2:'};
dlgtitle = 'ROI inputs:';
dims = [1 80];
definput = {'110','118'};
roi_input = inputdlg(prompt,dlgtitle,dims,definput)
a = str2double(roi_input);
a = a(1):a(2);
c = topcell(a,:)';
cell_location = cell((size(c)));
figure('Name','Calcium fluctuations of invivo cells'); clf;
 for i = 1:length(c)
     cell_location(:,i) = {['Cell ', int2str(c(:,i))]};
    plot(t , df_f(c(:,i),1:size(df_f,2)) + 1.5*(i-1)*ones(size(s(1,1:size(df_f,2)))),...
    'LineWidth',1.5);
    hold on; 
 
 end
ylabel('Calcium intensity');
xlabel('Time Frame');
legend(cell_location,'Location','best','NumColumns',2);
title(['Calcium fluctuations of invivo cells']);
hold off;
buttonname = 'no';
answer = 'no ';
while ~isequal(buttonname,'yes')
%include your switch-case code here

app = filterDesigner;
waitfor(app);
prompt = 'writh the variable name for SOS exported to workspace';
sos = input(prompt);
prompt = 'writhe the variable name for g scale exported to workspace';
g = input(prompt);
disp('coefficient is saved.');
y2 = []; 
for i = 1: size(df_f,1)
y2(i,:) = filtfilt(sos,G,df_f(i,:));
end 
%plotting| ROI 187 in M6 video. Finding the respected roi in denoised
%signal 
while ~isequal(answer,'no')
prompt = {' index number of ROI to display:'};
dlgtitle = 'ROI inputs:';
dims = [1 80];
definput = {'1'};
roi_input = inputdlg(prompt,dlgtitle,dims,definput)
kk = str2double(roi_input);
figure('Name',['df_f vs denoised cell' num2str(kk)]);clf; subplot(2,1,1); plot(df_f(kk,:));...
subplot(2,1,2); plot(y2(kk,:)); 
figure('Name','df_f vs denoised cell');clf; 
plot(df_f(kk,:))
hold on
plot(y2(kk,:),'LineWidth',2.5);
answer = questdlg('try another ROI?', ...
	'ROI', ...
	'yes','no','no');
% Handle response
switch answer
    case 'yes'
        disp([answer])
     
    case 'no'
        disp([answer])
     
end
end
buttonname = questdlg('accept the filter?', ...
	'filter', ...
	'yes','no','no');
% Handle response
switch buttonname
    case 'yes'
        disp([buttonname])
     
    case 'no'
        disp([buttonname])
     
end
end
%% Step 2 finding peaks and impulse train determination
%%findpeaks 
buttonname = 'no';
answer = 'yes';
while ~isequal(buttonname,'yes')
B_real = zeros(size(df_f)); 
%varry parameter_ to understand bette variation, please check findpeaks
%function from matlab 
%CHANGING THEM 

prompt = {'Threshold:','Min peak height:','Min distance between peaks:','Min prominence:','Min Width:','Max Width:'};
dlgtitle = 'FindPeaks inputs';
dims = [1 80];
definput = {'1e-5','0.5','5','0.5','3','50'};
findpeak_input = inputdlg(prompt,dlgtitle,dims,definput)
findpeak_input = str2double(findpeak_input);

%This will take about 3-5 min to process 
for i = 1: size(df_f,1)
        [pks, locs, width, p] = findpeaks(y2(i,:),'Threshold', findpeak_input(1), 'MinPeakHeight',...
        findpeak_input(2),'MinPeakDistance', findpeak_input(3),'MinPeakProminence',findpeak_input(4),'MinPeakWidth',findpeak_input(5),'MaxPeakWidth',findpeak_input(6));
        B_real(i,locs) = 1; 
  
end

% to investigate which ROIs spike and how many times? 
num_pks = sum(B_real.');
[rnk_pks,I] = sort(num_pks,'descend');
II = sort_cell(I); 
% Final array based on ROI number 
comb_rnks_pks = [rnk_pks ;II'];
% array based on ROIs that was above background
comb_rnks_2 = [rnk_pks; I]; 

%Check the most active cells 
while ~isequal(answer,'no')
prompt = {' index number of ROI to display:'};
dlgtitle = 'ROI inputs:';
dims = [1 80];
definput = {'1'};
roi_input = inputdlg(prompt,dlgtitle,dims,definput)
k = str2double(roi_input);
%Quality control based on visual plotting and checking movie
% at what time(s) spike
locs_new_peaks = find(B_real(k,:) ==1)
figure('Name',['Calcium fluctuations for cell ', num2str(k)]); clf; 
plot(y2(k,:));hold on;
plot(locs_new_peaks, y2(k,locs_new_peaks),'r*'); hold on;
ylabel('Calcium intensity');
xlabel('Time Frame');
title(['Calcium fluctuations for cell ', num2str(k)]);
legend('calcium signal', 'peaks'); hold off;
answer = questdlg('try another ROI?', ...
	'ROI', ...
	'yes','no','no');
% Handle response
switch answer
    case 'yes'
        disp([answer])
     
    case 'no'
        disp([answer])
     
end
end
% Different ways to show how active cells are
%number of peak for total cell:
zero_cell = size(find(rnk_pks ==0));
figure('Name','Frequency of real impulse'); clf; 
histogram(rnk_pks);
title('Frequency of real impulse'); 

[row,col] = find(B_real ==1);
[wcell,ind] = sort(row,'ascend');
%time of spiking from all active ROIs
tim_spk = col(ind);
% time of spiking + number of ROI that is involved
comb_cor = [sort_cell(wcell) tim_spk];
% histogram of time of spiking
figure('Name', 'histogram of time of spiking'); clf; h2= histogram(tim_spk, 'BinWidth',30);
title(' histogram of time of spiking');
% number of spikes + which active ROI ( only for active cells) 
cell_spike = comb_rnks_2(:,find(rnk_pks~=0)); 
buttonname = questdlg('Accept peaks?', ...
	'filter', ...
	'yes','no','no');
% Handle response
switch buttonname
    case 'yes'
        disp([buttonname])
     
    case 'no'
        disp([buttonname])
     
end
end
%% Step 3
% generating generated cells
buttonname = 'no';
while ~isequal(buttonname,'yes')
 disp('Choose rnk_pks variable in distributionFitter app.')
app = distributionFitter; 
disp('close the distributionFitterapp once youre done');
waitfor(app);
prompt = 'type the name that saved fitted distribution';
pd_np = input(prompt)
% One-sample Kolmogorov-Smirnov test/ a nonparametric test
%'CDF' — cdf of hypothesized continuous distribution
% if f is 1 this indicates the rejection of the null hypothesis at the Alpha significance level.
% p shows the p value 
%cv ~ critical value 
[f, p, ~,~] = kstest(rnk_pks.', 'CDF',pd_np)

disp('if f is 1 this indicates the rejection of the null hypothesis at the Alpha significance level.')
disp('p shows the p value ')
othergrid_rnk = transpose(linspace(min(rnk_pks(:)),max(rnk_pks(:)),100));
figure('Name','Histogram of Number of impulses for each ROIs in real data with pdf');clf;
histogram(rnk_pks,25,'Normalization','pdf','FaceColor','g'); %plot original data
w_rnk = pdf(pd_np,othergrid_rnk);
hold on
plot(othergrid_rnk,w_rnk,'LineWidth',2,'Color','r') %plot GMM over original data
ylabel('Frequency'); 
xlabel('number of peaks'); 
title('Histogram of Number of peaks for each ROIs');
legend('Real data,number of peaks');
hold off
buttonname = questdlg('Accept probability distribution?', ...
	'filter', ...
	'yes','no','no');
% Handle response
switch buttonname
    case 'yes'
        disp([buttonname])
     
    case 'no'
        disp([buttonname])
     
end
end
% truncate probability dist object to the number of max impulse train sum which is 8 for this case 
pd_np= truncate(pd_np, 0,max(rnk_pks));
%random number of peaks based on negative binomial fit with same number of
%active ROI - as an example 
rnd_np = random(pd_np,size(df_f,1),1);
[nreal, edgesreal] = histcounts(rnk_pks,8);
%1000 random number of peaks based on negative binomial fit with same number of
%active ROI 
rnd_imp = zeros(size(rnk_pks,2),1000); 
n = 0 ; 
if n < size(rnk_pks,2) 
for i =1 :1000
    j = 0;
j = random(pd_np,size(rnk_pks,2),1);
rnd_imp(:,i) = j; 

 
end
    [n, edges,bin] = histcounts(rnd_imp);
    n = round(n./1000);
end
if sum(n) > size(df_f,1)
    a = abs(size(df_f,1) - sum(n));
    n(:,1)= n(:,1) -a; 
    else 
      a = abs(size(df_f,1) - sum(n));
    n(:,1)= n(:,1) + a;
end
% n is number of cells within each spikes from 0 .....8 
%bin_spks an array of numbe of spikes with an appropriate of  number of cells 
bin_spks= [0:size(n,2)-1;n];
figure('Name','Frequency of B real and B generated'); clf;
subplot(2,1,1); histogram(rnk_pks, max(rnk_pks),'FaceColor','g'); 
title('real number of peaks within each active ROIs');
subplot(2,1,2); histogram(rnd_np,max(rnk_pks),'FaceColor',[0.75 0.75 0.75]);
title('random number of peaks within each ROIs NegativeBinomial dist');

%Giving cells random number of spikes based on probability of distribution
 art_spk =[];
  h= 0; 
   for j = 1:size(bin_spks,2) 
      
       for m = h+1:h+bin_spks(2,j) 
          art_spk(m,:) = bin_spks(1,j);  
       end
       h = h + bin_spks(2,j);
   end
%Random permutation 
rng('shuffle');
rand_perm_spk = randperm(size(art_spk,1));
art_spk2 = art_spk; 
art_spk2(rand_perm_spk,:) = art_spk(:,:);

%%Generate a new generated binary matrix 
k = 0;
% this will take about 2-3 min since it's generating 100 binary matrix 
B_gen = zeros(size(rnk_pks,2), length(B_real),100);
for j = 1:100
for i = 1:(size(rnk_pks,2))
    k= 0 ;
    if (art_spk2(i,:) ~=0)
        k = art_spk2(i);
       rng('shuffle')
       rnd_numpks = randi([1 600],1,k);
        B_gen(i,rnd_numpks,j) = 1;
    end
end 
end 

% where the peak happened? 
[~,wpks] = find(B_real==1);
[a,bb,b] = histcounts(wpks,20);
[~,wpks_art] = find(B_gen==1);
[n, edges] = histcounts(wpks_art,20);
n = round(n./100); 
edges = round(edges./100);
%showing a  histogram of both real data(green) and 100 generated
%data(grey) 
figure('Name','Frequency of impulses in time for both real and generated data');clf; histogram('BinCounts',n,'BinEdges',edges,'FaceColor', [0.75 0.75 0.75]); hold on ;
histogram('BinCounts',a,'BinEdges',bb,'FaceColor','g');
%%defining threshold from generated data
%generating threshold for counting S based on generated impulses 

%% Output
%This step requires an intervention of users. 1. investigating the optimum
%window size, 2. Paying attention to optimal K number and plot based on
%that. 
%%INPUT: 
 %optimum w
 % finding maximum active number within a size of window /w 
 buttonname = questdlg('Choose optimal window size?', ...
	'filter', ...
	'yes','no ','no ');
% Handle response
switch buttonname
    case 'yes'
        disp([buttonname])
     
    case 'no'
        disp([buttonname])
     
end

if buttonname  =='yes'
Sth = [];
prompt = {'z score:'};
dlgtitle = 'z score input';
dims = [1 30];
definput = {'1.28'};
z_score = inputdlg(prompt,dlgtitle,dims,definput);
z = str2double(z_score)
w_sz= size(df_f,2)/6;
for sz_bin= 1:w_sz
     t_art = []; 
 for j3 = 1:100
for j =  1:(size(df_f,2)-sz_bin);
    [locs_pks,col] = find(B_gen(:,j:j+sz_bin,j3) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
     t_art(:,j,j3) = m ;   
     
end
 end
  [N2, edges2] = histcounts(t_art);
N2 = round(N2./100);
center_art = edges2 - 0.5;
 center_art(:,1) = []; 
comb_art = [center_art;N2];
 a_art = find(comb_art(2,:)~=0);
Sth(sz_bin) = (z*mean(comb_art(1,a_art))/sqrt(length(comb_art(1,a_art))) + mean(comb_art(1,a_art)));
end
 act_p = [];
for sz_bin = 1:w_sz
 t_bin = []; 

for  j = 1:(size(y2,2) - sz_bin);
    [Cell, ~] = find(B_real(:,(j:j+sz_bin)) ==1); 
    unique_cell = unique(Cell);
    m = size(unique_cell,1);
    t_bin(:,j) = m ; 
end
  act_p(sz_bin,:) = (max(t_bin)- Sth(sz_bin))/sz_bin;
end
opt_w = [];
while isempty(opt_w)
prompt = {'Threshold:','Min peak prominence:'};
dlgtitle = 'peak input for optimized peak';
dims = [1 50];
definput = {'1e-5','0.5'};
peak_input = inputdlg(prompt,dlgtitle,dims,definput);
peak_input = str2double(peak_input)
 [p,l, w,pr] = findpeaks(act_p, 'Threshold', peak_input(1),'MinPeakprominence',peak_input(2));
 figure('Name', 'Opt window size'); clf; 
plot(act_p);hold on;
plot(l, act_p(l,:),'r*'); hold on;
% set(gca,'ytick',[])
ylabel('Max(active cell /w');
xlabel('Window size ');
legend('calcium signal', 'peaks'); hold off;
title('optimum window size');
[~,i] = max(l./p);
opt_w = l(i)
 opt_w2 = opt_w
end
else 
 prompt = {'manual w:'};
dlgtitle = 'manual window input';
dims = [1 30];
definput = {'6'};
manual_input = inputdlg(prompt,dlgtitle,dims,definput);
opt_w2 = str2double(manual_input)
end 
 %%counting S for both real and generated  with a fixed w of 5 time frame
 %%(sample)
% counting S in generated binary matrix with fixed window size
prompt = {'beginning time:'};
dlgtitle = 'beginning time (tau) input';
dims = [1 60];
definput = {'1'};
time_input = inputdlg(prompt,dlgtitle,dims,definput);
tau = str2double(time_input)

  s_gen = []; 
 for j3 = 1:size(df_f,2)/6
for j =  tau:(size(y2,2)-opt_w2);
    [locs_pks,col] = find(B_gen(:,j:j+opt_w2,j3) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
     s_gen(:,j,j3) = m ;   
end
 end
   
% counting S in real binary matrix with fixed window size
 
s_real = []; 
 
for  j = tau:(size(y2,2) - opt_w2);
    [locs_pks, col] = find(B_real(:,(j:j+opt_w2)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    s_real(:,j) = m ; 
          
end
 
 [N2, edges2] = histcounts(s_gen);
N2 = round(N2./100);
[N, edges] = histcounts(s_real);
center_art = edges2 - 0.5;
 center_art(:,1) = []; 
comb_art = [center_art;N2];
center_real =edges - 0.5; 
center_real(:,1) = [];
comb_real = [center_real;N]; 
 a_real= find(comb_real(1,:)~=0);
 a2_real = max(comb_real(1,a_real));
   a_art = find(comb_art(2,:)~=0);
   a2_art = max(comb_art(1,a_art));
   % threshold using generated cells 
   if buttonname == 'yes'
   figure('Name','S real and Threshold Sth');clf; bar(s_real); hold on; plot(1:length(s_real), ones(1,length(s_real))*Sth(opt_w));
 title('S real and Threshold Sth'); 
 delta_s = (s_real - Sth(opt_w2)) ./ opt_w2; 
delta_s(delta_s <0 )= 0; 
   else
      Sth = [];
prompt = {'z scpore:'};
dlgtitle = 'z score input';
dims = [1 30];
definput = {'1.28'};
z_score = inputdlg(prompt,dlgtitle,dims,definput);
z = str2double(z_score)
       Sth = (z*mean(comb_art(1,a_art))/sqrt(length(comb_art(1,a_art))) + mean(comb_art(1,a_art)))
       delta_s = (s_real - Sth); 
delta_s(delta_s <0 )= 0; 
          figure('Name','S real and Threshold Sth');clf; bar(s_real); hold on; plot(1:length(s_real), ones(1,length(s_real))*Sth);
 title('S real and Threshold Sth'); 
   end 
%countS with a fixed threshold

 figure('Name','distribution of synchroncity'); clf; 
 h3= histogram('BinCounts',N,'BinEdges',edges,'FaceColor','g');hold on;

 title(['distribution of synchroncity, size of bin=',int2str(opt_w2)]); hold on; 

 h4= histogram('BinCounts',N2,'BinEdges',edges2,'FaceColor',[0.75 0.75 0.75]);

legend('Real data','Null data','location','northeast');
ylabel(['#of ROIs in time windoe of',int2str(69)]); 
xlabel(['number of spikes for all ROIs within time bin of ',int2str(69)]);
hold off;

[~,time] = find(delta_s~=0);
j = 0; 
t_imp= [];
for i = 1: length(time)-1
    if time(i+1) - time(i) >4*opt_w2
    j = j +1;
    t_imp(j) = time(i);
    end
end
num_com = j;
for j = 1: num_com
%%cross correlation between space(X and Y coordinate) and time of
twindow = [];
if j ==1 
twindow(:,1) = time(1): t_imp(j)  ; 
twindow(:,2) = time(1) + opt_w2 : t_imp(j) + opt_w2; 

else 
    clear twindow 
  twindow(:,1) = t_imp(j-1): t_imp(j)  ; 
twindow(:,2) = t_imp(j-1) + opt_w2 : t_imp(j) + opt_w2;   
end
% which ROIs has spike within this range 
sync_roi =  cell(size(twindow,1),1);
for i = 1: size(twindow,1)
    t = twindow(i,:); 
    [roi, ~] = find(B_real(:,t(1):t(2))==1); 
    sync_roi{i,1} = roi.'; 
end

sync_roi_vct = horzcat(sync_roi{:});
[sync_roi_vct,~,~] = unique(sync_roi_vct,'first');
sync_roi_vct2 = sort(sync_roi_vct);
sync_roi_vct = sort(sort_cell(sync_roi_vct)); 
sync_roi_str = string(sync_roi_vct);
sync_roi_loc= s_loc(:,sync_roi_vct.');

figure('Name','centroid of synchronous cells' ) ;clf; scatter(sync_roi_loc(1,:), sync_roi_loc(2,:),'green','filled');
dx = 0.3; dy = 0.3; % displacement so the text does not overlay the data points
text(sync_roi_loc(1,:)+dx,sync_roi_loc(2,:)+dy, sync_roi_str);
% xlim([ 0 636.16]);
% ylim([0 636.16]);
title 'centroid of synchronous cells' 


        rng('default') % For reproducibility
        eva = evalclusters( sync_roi_loc.','kmeans','silhouette','klist',[1:12])
figure('Name','cluster evaluation'); clf;
plot(eva);
%Unsupervised k-mean clustering based on evaluation(
%k = 4) 

opts = statset('Display','final');
[idx,C,sumd,D] = kmeans(sync_roi_loc.',eva.OptimalK,'Distance','sqeuclidean',...
    'Replicates',5);
a = sync_roi_loc.';

% out put data : 1. time of communication 2. a matrix of centroids
% information of synchrounous cell with its designated id number of
% cluster. the size is x*3 which x shows the number of synchoronous cells
%%desire directory to save figures and outputs 
disp('Choose/create a output folder')
output_path = uigetdir
output_path = append(output_path,'/');
csvwrite(append(output_path,'time of communication_',string(j),'.csv'),twindow) ;
q_met = [];
for i = 1: eva.OptimalK
    if any(sumd(i)) == 0
        q_met(i,:) = 0; 
    else
q_met(i,:)= length(a(idx==i,1))/sumd(i);
    end
    end 
[~,syn_c] = max(q_met)

figure('Name','quality metric');clf; plot(q_met); hold on; 
plot(syn_c, q_met(syn_c),'r*'); hold on; 
ylabel('quality metric for cells that are communicating');
xlabel('number of cluster');

 cell_comm = idx ==syn_c;
 cell_comm = a .* cell_comm;
 cell_comm = cell_comm(cell_comm ~= 0); 
 cell_comm = reshape(cell_comm, length(cell_comm)/2,2);
final_feature = array2table([cell_comm]);
final_feature.Properties.VariableNames = {'CentroiX', 'CentroidY'};
writetable(final_feature,append(output_path,'centroid_synchrouns_cell_',string(j),'.csv'));
end

%saving figures as jpg 
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');

  saveas(FigHandle, fullfile(output_path, [FigName, '.jpg']));
end