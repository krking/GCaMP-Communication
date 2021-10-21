 %% defining threshold from generated data -Figure 5j
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


th= [ min(min(min(min(count_s3(:,:,:,:))))) :0.05: max(max(max(max(count_s3(:,:,:,:)))))];
% Figure 5J
%number of communication based on window size
th_cl = []; 

for sz_bin =2:15 
b = (Sth(sz_bin)*mean(active_roi_2))/size(m_cell_general,1);

    m=permute(count_s3(:,:,sz_bin,:), [1 2 4 3]);
% for ss= 1:length(th)
k = 1;

%   m = flip(m,3);
%  np3 = padarray(m, [1 1],0,'both');

%you can vary the percentage- look at 5J
 th = 0.75* max(max(max(count_s3(:,:,sz_bin,:))));

n=(m>th); %make binary if above threshold

np = padarray(n,[1 1],0,'both'); %pad w 0's
 np2 = padarray(n, [1 1],0,'both');

np=double(np);

sn=size(n); % size of original binary matrix

snp=size(np); % size of padded binary matrix

 

kern = [1,1,1;1,0,1;1,1,1];

kern_cmp = bitcmp(kern);

c = 1; % cluster index   

box_com2 =[];
for t = 1:sn(3)

    % must first expand locations that were pushed forward and are >1

    [rowList,colList] = find(n(:,:,t)>1); 

    ds1 = size(rowList);

    while (length(rowList)>0)

        [p,row,col]=propagate(np((rows(1)-1):(rows(1)+1),(cols(1)-1):(cols(1)+1),t),i,j);

        ds2 = size(row);

        np((rowList(1)-1):(rowList(1)+1),(colList(1)-1):(colList(1)+1),t) = p;

        if rowList>0

            rowList = [rowList(2:end); row]; % remove first element, add new locs

            colList = [colList(2:end); col]; % remove first element, add new locs

            %ds3 = size(rowList)

            %ds4 = size(colList)

        end

    end

    

    for i=2:(sn(1)+1) 

        for j=2:(sn(1)+1)

            % create a mini 3x3 

            x3mtx = np((i-1):(i+1),(j-1):(j+1),t);

            % skip this loc if np(i,j,t) is 0 or has already been set >1;

            if np(i,j,t)==1 

                % Option 1: convert middle to max of surrounding nonzeros

                a = max(max(x3mtx((kern.*x3mtx)>1))); % max of surrounding >1 values

                if a > 1

                    % set middle to max of outside values >1

                    np(i,j,t)=a;

                   [t,i,j,np(i,j,t)]            

                % Option 2: set clusterNum, increment clusterCounter,

                % then propagate

                else

                    np(i,j,t) = c;

                    rowList=i;

                    colList=j;

                    while (length(rowList)>0)

                        [p,row,col]=propagate(np((i-1):(i+1),(j-1):(j+1),t),i,j);

                        np((rowList(1)-1):(rowList(1)+1),(colList(1)-1):(colList(1)+1),t) = p;

                        if length(row>0)

                            rowList = [rowList; row];

                            colList = [colList; col];

                        end

                        ds5 = rowList;

                        ds6 = row;

                        ds7 = colList;

                        ds8 = col;

                        % remove first element of rowList and colList

                        if length(rowList==1)

                            rowList=[];

                            colList=[];

                        else

                            rowList = rowList(2:end); 

                            colList = colList(2:end); 

                        end

                        %pause

                    end

                    np(i,j,t);

                    c=c+1 % increment cluster num

                  box_com2(c-1,:)= [t,i,j,np(i,j,t)];

%                     figure(2);
% 
%                     imagesc(np(:,:,t));

                    %pause 

                end

            end

        end

    end

%     figure(1);
% 
%     imagesc(np(:,:,t));
% 
%     title("Time " + t + " - Cluster " + c)
% 
%     c_mov(k) = getframe(gcf);

    %k=k+1;

    % Push >1 values forward in time

    if t<sn(3)

        %d1 = np(:,:,t)

        %d2 = np(:,:,(t+1))

        np(:,:,(t+1)) = max(np(:,:,t),np(:,:,(t+1))).*np(:,:,(t+1));

        %d3 = np(:,:,t)

        %d4 = np(:,:,(t+1))

    end

end
 b = unique(np(np2==1));
% th_cl(sz_bin,:) = [th max(b)];
%   end

%   a=  th_cl(2:end,2);
% a = a - 0.4; 

m = 1;

if size(box_com2,1)>1
for i  =1:size(box_com2,1)-1
       if (box_com2(i,2)~=box_com2(i+1,2)) && (box_com2(i,3)~=box_com2(i+1,3))
              
     m = m+1;
       end
end
end

th_cl(sz_bin) = m;
end
 