clear all 
close all

a=rand(10,1)
l1_norm_a=norm(a, 1);
l2_norm_a=norm(a, 2);
s=a;
k1=0.5*l1_norm_a;
k2=l2_norm_a;

%try k1 such that it acheive sparsenesss=0.6
k1=k2*(sqrt(numel(a))-0.6*(sqrt(numel(a))-1));

sparseness=(sqrt(numel(a))-k1/k2)/(sqrt(numel(a))-1)

nn=1;
[v,usediters] = projfunc( s, k1, k2, nn )
sparseness=(sqrt(numel(v))-norm(v,1)/norm(v,2))/(sqrt(numel(v))-1)

ind=1;

for i=1:0.1:6  
    
    two_col_mat(ind, 1)=i;
    
    two_col_mat(ind, 2)=(sqrt(10)-i)/(sqrt(10)-1)
    
    ind=ind+1;
    
end

