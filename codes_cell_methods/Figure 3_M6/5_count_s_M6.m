%% Method for quantifying normalized number of spikes Fig 3E and G
%% defining threshold from generated data -Figure 3F&G
% depends on how many window size you wanna investigate the for loops changes
% this will take about 10-15 min 
Sth = [];
w_sz= 100;
for sz_bin= 1:w_sz
     t_art = []; 
 for j3 = 1:100
for j =  1:(size(df_f,2)-sz_bin);
    [locs_pks,col] = find(art_cell(:,j:j+sz_bin,j3) ==1); 
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
Sth(sz_bin) = (1.96*mean(comb_art(1,a_art))/sqrt(length(comb_art(1,a_art))) + mean(comb_art(1,a_art)));
end



% counting S 
count_s = [];
for sz_bin = 1:w_sz
    
 t_bin = []; 

 
for  j = 1:(size(y2,2) - sz_bin);
    [locs_pks, ~] = find(m_cell(:,(j:j+sz_bin)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    t_bin(j+round(sz_bin/2)) = m;  

    if  t_bin(j) > Sth(sz_bin) 
    count_s(sz_bin,j) = (t_bin(j) - (Sth(sz_bin)))./sz_bin;
else 
    count_s(sz_bin,j) = NaN; 
end

end

end
count_s(count_s ==0) = NaN;


figure(1);clf;

[nr,nc] = size(count_s);
 h = pcolor([count_s nan(nr,1); nan(1,nc+1)]);
 colormap(jet); 
 shading flat;
 set(gca, 'ydir', 'reverse');
 
 % finding optimum tau 
 
 [i,ii] = max(count_s,[],2); 
a = isnan(i);
i (isnan(i)) = [];
ii ( a==1) =[];
 
  j = i > (1.282*mean(i)/sqrt(length(i)) + mean(i)); 
 i = i( j==1);
 ii = ii(j ==1);
 ii = sort(ii);
 [i3,ia,~] = unique(ii); 
 %optimum tau 
 t_consider = i3; 
 
 %optimum window 
 count_s_imp= mean(count_s(:,t_consider(1):t_consider(end)),2);
count_s_imp= mean(count_s_imp,2);

 [~,opt_w] = max(count_s_imp,[],1)
 % optimum size for window is 12 
 %% counting S for both real and artificial -Figure 3E
% counting S in artificial binary matrix with fixed window size
  t_art = []; 
 for j3 = 1:100
for j =  1:(size(y2,2)-opt_w);
    [locs_pks,col] = find(art_cell(:,j:j+opt_w,j3) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
     t_art(:,j,j3) = m ;   
end
 end
   
% counting S in real binary matrix with fixed window size
 
 t_bin = []; 
 
for  j = 1:(size(y2,2) - opt_w);
    [locs_pks, col] = find(m_cell(:,(j:j+opt_w)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    t_bin(:,j) = m ; 
          
end


 
 [N2, edges2] = histcounts(t_art);
N2 = round(N2./100);
[N, edges] = histcounts(t_bin);
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
 figure(2);clf; bar(t_bin); hold on; plot(1:length(t_bin), ones(1,length(t_bin))*Sth(opt_w));
%countS with a fixed threshold

 figure(3); clf; 
 h3= histogram('BinCounts',N,'BinEdges',edges,'FaceColor','g');hold on;

 title(['distribution of synchroncity, size of bin=',int2str(opt_w)]); hold on; 

 h4= histogram('BinCounts',N2,'BinEdges',edges2,'FaceColor',[0.75 0.75 0.75]);

legend('Real data','Null data','location','northeast');
ylabel(['#of ROIs in time windoe of',int2str(69)]); 
xlabel(['number of spikes for all ROIs within time bin of ',int2str(69)]);
hold off;
