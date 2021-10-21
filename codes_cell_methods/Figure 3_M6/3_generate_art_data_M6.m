figure(1); clf; 
histogram(rnk_pks)
%Fit probability distribution object to experimental data for number of
%peaks using negative binomial 
pd_np = fitdist(rnk_pks.', 'NegativeBinomial');

% One-sample Kolmogorov-Smirnov test/ a nonparametric test
%'CDF' — cdf of hypothesized continuous distribution
% if f is 1 this indicates the rejection of the null hypothesis at the Alpha significance level.
% p shows the p value 
%cv ~ critical value 
[f, p, ss,cv] = kstest(rnk_pks.', 'CDF',pd_np)
othergrid_rnk = transpose(linspace(min(rnk_pks(:)),max(rnk_pks(:)),100));
figure(1);clf;
histogram(rnk_pks,25,'Normalization','pdf','FaceColor','g'); %plot original data
w_rnk = pdf(pd_np,othergrid_rnk);
hold on
plot(othergrid_rnk,w_rnk,'LineWidth',2,'Color','r') %plot GMM over original data
ylabel('Frequency'); 
xlabel('number of peaks'); 
title('Histogram of Number of peaks for each ROIs');
legend('Real data,number of peaks','NegativeBinomial fit');
hold off
% truncate function 
pd_np= truncate(pd_np, 0,8);
%random number of peaks based on negative binomial fit with same number of
%active ROI - as an example 
rnd_np = random(pd_np,size(df_f,1),1);
[nreal, edgesreal] = histcounts(rnk_pks,8);
figure(2); clf;
subplot(2,1,1); histogram(rnk_pks, 8,'FaceColor','g'); 
title('real number of peaks within each active ROIs');
subplot(2,1,2); histogram(rnd_np,8,'FaceColor',[0.75 0.75 0.75]);
title('random number of peaks within each ROIs NegativeBinomial dist');
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
end
% n is number of cells within each spikes from 0 .....8 
%bin_spks an array of numbe of spikes with an appropriate of  number of cells 
bin_spks= [0:size(n,2)-1;n];
figure(3); clf;
subplot(2,1,1); histogram(rnk_pks, 'BinEdges',edges,'FaceColor','g'); 
title('real number of peaks within each ROIs');
subplot(2,1,2); h = histogram('BinEdges',edges,'BinCounts',n,'FaceColor',[0.75 0.75 0.75]);
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
%% Generate a new artificial binary matrix 
k = 0;
% this will take about 2-5 min since it's generating 100 binary matrix
art_cell = zeros(size(rnk_pks,2), length(m_cell),100);
for j = 1:100
for i = 1:(size(rnk_pks,2))
    k= 0 ;
    if (art_spk2(i,:) ~=0)
        k = art_spk2(i);
       rng('shuffle')
       rnd_numpks = randi([1 600],1,k);
        art_cell(i,rnd_numpks,j) = 1;
    end
end 
end 

% where the peak happened? 
[~,wpks] = find(m_cell==1);
[a,bb,b] = histcounts(wpks,20);
[~,wpks_art] = find(art_cell==1);
[n, edges] = histcounts(wpks_art,20);
n = round(n./100); 
edges = round(edges./100);
%showing a  histogram of both real data(green) and 100 artificial
%data(grey) 
figure(4);clf; histogram('BinCounts',n,'BinEdges',edges,'FaceColor', [0.75 0.75 0.75]); hold on ;
histogram('BinCounts',a,'BinEdges',bb,'FaceColor','g');
% Set up a threshold(red) based on average bin count of 100 artificial data
figure(5); clf;  histogram('BinCounts',a,'BinEdges',bb,'FaceColor','g'); hold on; ...
    plot(1:600, ones(1,600)*(mean(n)+2*std(n))); hold off;
