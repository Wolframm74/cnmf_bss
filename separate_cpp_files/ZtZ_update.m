function [Z_ol]=ZtZ_update(Z_ol)

display('inside ZtZ_update');

[O_static, L_static]=size(Z_ol);

R_z=transpose(Z_ol)*Z_ol;

% R_z_indices=find(R_z>1*10^(-9));
% 
% R_z(R_z_indices)=1*10^(-9);



R_z_mask=zeros(L_static, L_static);

R_z_mask=diag(ones(L_static, 1), 0);

for i=1:(L_static-2)

R_z_mask=R_z_mask+diag(exp(-i*2)*ones(L_static-i, 1), i);

R_z_mask=R_z_mask+diag(exp(-i*2)*ones(L_static-i, 1), -i);

end

        figure, 
        imagesc(Z_ol);
        title('MEX: Z_ol: before');

        figure,
        imagesc(R_z);        
        title('MEX: ZtZ_update: before');
        
        figure,
        imagesc(R_z_mask);
        title('MEX: R_z_mask');        
        
Z_ol=ZtZ_update_local(Z_ol, R_z.*R_z_mask, O_static, L_static);        

R_z=transpose(Z_ol)*Z_ol;

        figure, 
        imagesc(Z_ol);
        title('MEX: Z_ol: after');

        figure,
        imagesc(R_z);        
        title('MEX: ZtZ_update: after');
 

end

function [Z_ol]=ZtZ_update_local(Z_ol, R_z, O_static, L_static)

[R_z]=ZtZ_update_sort_diagonal( R_z,  L_static);

gamma=2;
%gamma=30;

figure, 
imagesc(R_z);
title('MEX: inside ZtZ_update_local: R_z.*R_z_mask');

onesmat_OL=ones(O_static, L_static);

%Z_ol=Z_ol.*((1-(1/gamma))*onesmat_OL+(1/(2*gamma))*(Z_ol*R_z+Z_ol*transpose(R_z))./(Z_ol*transpose(Z_ol)*Z_ol));
Z_ol=Z_ol.*((1-(1/gamma))*onesmat_OL+(1/(2*gamma))*(Z_ol*R_z+Z_ol*transpose(R_z))./(2*Z_ol*transpose(Z_ol)*Z_ol));


end

%find the indices of the L_target=3 most dominant diagonal elements.
%compute their average amplitude/energy
%Set R_z at these indices equal to the average: a crude normalization. 
function [R_z]=ZtZ_update_sort_diagonal( R_z,  L_static)

L_target=3;

diagvec=diag(R_z);

Ltargetx2_info_colmatrix=zeros(L_target, 2);

%find the indices of the L_target=3 most dominant diagonal elements.
for l=1:L_target

[amp, index]=max(diagvec);    
    
Ltargetx2_info_colmatrix(l, 1)=amp;

Ltargetx2_info_colmatrix(l, 2)=index(1);

%zero it but not to worry it will be reset to nonzero at some point
diagvec(index(1))=0;
    
end

%sum their amplitudes
total_energy=sum(squeeze(Ltargetx2_info_colmatrix(:, 1)));

%compute their average amplitude/energy
avg_energy=total_energy/L_target;

%Halve all energies in R_z. Over time this should persist only the
%L_target=3 most significant clusters and attenuate the rest...
%R_z=0.5*R_z;

max_diagvec=Ltargetx2_info_colmatrix(1,1);

diagvec=shaping_helper_1(diagvec, max_diagvec, L_static);

%use the 2nd column of the info column matrix to index and set the
%avg_energy value. 
diagvec(squeeze(Ltargetx2_info_colmatrix(:, 2)))=avg_energy;

%Set R_z at these indices equal to the average: a crude normalization. 
R_z(logical(eye(L_static))) = diagvec;

end

function [diagvec]=shaping_helper_1(diagvec, max_diagvec, L_static)

elementwise_scaling_vec=zeros(L_static, 1);

elementwise_scaling_vec=ones(L_static,1)-exp(-(10/(2*max_diagvec))*diagvec);

diagvec=diagvec.*elementwise_scaling_vec;

end