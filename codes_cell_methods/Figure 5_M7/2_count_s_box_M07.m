%% IF your movie has a high synchrony you need split movies

%max diameter of ROI = 9 um 

%collecting rois based on box
box = strcat('box');
a3 =  0:(round(max(s_loc(1,:))-50))/25;
a4 = [25*a3 ; 50+25*a3]; 
empt_string = strings(1,length(a3));
j3= [];
for i = 1:length(a3)
    for j = 1:length(a3)
         j3{i,j} = struct;
        [~,roi] = find(s_loc(1,:) <50+(25*a3(j)) & s_loc (1,:)>25*a3(j));
        loc_new = s_loc(:,roi);
        [~,new_roi] = find(loc_new(2,:)< 50+(25*a3(i)) & loc_new(2,:)>25*a3(i));
             empt_string(:,i) =  strcat('box', string(i),'box',string(j));
     j3{i,j}.(empt_string(:,i)) = roi(new_roi);
    end
end

total_roi = []; 
    active_roi =[];
for i = 1: length(j5)
    for ii =1 : length(j5)
     roi = struct2array(j5{i,ii});
        [roi2,~] = find(sort_cell ==roi); 
        m_cell_sub = m_cell(roi2,:);
        [roi3,~] = find(m_cell_sub ==1);
        roi3 = unique(roi3);
         roi3 = unique(roi3);
active_roi(i,ii) = length(roi3);
 total_roi(i,ii) = length(roi);
        
    end 
end

figure(2);clf; subplot(1,2,1); histogram(active_roi); subplot(1,2,2); histogram(total_roi);
active_roi_2 = reshape(active_roi,[],1);
active_roi_2(isnan(active_roi_2)) = [];
total_roi_2 = reshape(total_roi,[],1);

j5 = j3; 

count_s3 = []; 
for sz_bin= 1:w_sz
for i = 1: length(j5)
    for ii =1 : length(j5)
     roi = struct2array(j5{i,ii});
        [roi2,~] = find(sort_cell ==roi); 
        m_cell_sub = m_cell(roi2,:);
        [roi3,~] =find(m_cell_sub==1);
        roi3 = unique(roi3);
thr_real = (Sth(sz_bin)*mean(active_roi_2))/size(m_cell_general,1);

 t_bin = [];  
for  j = 1:(size(y2,2) - sz_bin);
    [locs_pks, ~] = find(m_cell_sub(:,(j:j+sz_bin)) ==1); 
    unique_cell = unique(locs_pks);
    m = size(unique_cell,1);
    t_bin(j+round(sz_bin/2)) = m; 
    
count_s3(i,ii,sz_bin,j) = (t_bin(j) - thr_real)./sz_bin;

end
    end

end
end
count_s3( count_s3<=0) =NaN;
 count_s4 = count_s3; 
 
  %Figure 5-D
  figure(1);clf;
% size of 4D array is a box of 11*11 with a 100 varied window size and time
% scale 
  a = permute(count_s3(:,:,5,458),[1 2 3 4]);
  [nr,nc] = size(a);
 h = pcolor([a nan(nr,1); nan(1,nc+1)]);
 colormap(parula); 
 shading flat;
 caxis([0.1 1]);
 %Figure 5-E
 figure(2);clf;
  a = permute(count_s3(:,:,5,462),[1 2 3 4]);
  [nr,nc] = size(a);
 h = pcolor([a nan(nr,1); nan(1,nc+1)]);
 colormap(parula); 
 shading flat;
 caxis([0.1 0.7]);
 