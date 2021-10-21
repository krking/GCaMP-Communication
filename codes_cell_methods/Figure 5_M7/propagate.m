function [pm, r, c] = propagate(m,i,j)

 

% m= 3x3 matrix to be propagated to boxes that equal 1 

% pm= propagated 3x3 matrix

% add r & c as to list of centroids to be propagated

pm=m;

vmax = max(max(m));

[r,c]=find(m==1);

pm(m==1)=vmax;

r=r+i;

c=c+j;