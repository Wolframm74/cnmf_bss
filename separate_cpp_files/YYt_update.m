function [Y_lk]=YYt_update(Y_lk)

[L_static, K_static]=size(Y_lk);

R_y=Y_lk*transpose(Y_lk);

% R_y_indices=find(R_y>(1000));
% 
% R_y(R_y_indices)=(1000);

R_y_mask=zeros(L_static, L_static);

R_y_mask=diag(ones(L_static, 1), 0);

for i=1:(L_static-2)

R_y_mask=R_y_mask+diag(exp(-i*2)*ones(L_static-i, 1), i);

R_y_mask=R_y_mask+diag(exp(-i*2)*ones(L_static-i, 1), -i);

end

%         figure,
%         imagesc(R_y);
%         title('MEX: YYt_update');
        
        figure, 
        imagesc(Y_lk);
        title('MEX: Y_lk: before');
        
        figure,
        imagesc(R_y);
        title('MEX: YYt_update: before');
        
        figure,
        imagesc(R_y_mask);
        title('MEX: R_y_mask');

Y_lk=YYt_update_local(Y_lk, R_y.*R_y_mask, L_static, K_static);        
R_y=Y_lk*transpose(Y_lk);

        figure, 
        imagesc(Y_lk);
        title('MEX: Y_lk: after');

        figure,
        imagesc(R_y);
        title('MEX: YYt_update: after');
        
        %V_nk=transpose(Y_lk);
        
end


function [Y_lk]=YYt_update_local(Y_lk, R_y, L_static, K_static)

gamma=4;
%gamma=30;

[R_y]=YYt_update_sort_diagonal( R_y,  L_static);

        figure,
        imagesc(R_y);
        title('MEX: inside YYt_update_local: R_y.*R_y_mask');

onesmat_LK=ones(L_static,K_static);

%Y_lk=Y_lk.*((1-(1/gamma))*onesmat_LK+(1/(2*gamma))*(transpose(R_y)*Y_lk+R_y*Y_lk)./(Y_lk*transpose(Y_lk)*Y_lk));
Y_lk=Y_lk.*((1-(1/gamma))*onesmat_LK+(1/(2*gamma))*(transpose(R_y)*Y_lk+R_y*Y_lk)./(2*Y_lk*transpose(Y_lk)*Y_lk));

end

function [Y_lk]=YYt_update_local_EG(Y_lk, R_y, L_static, K_static)

aeta=0.1;

[R_y]=YYt_update_sort_diagonal( R_y,  L_static);

        figure,
        imagesc(R_y);
        title('MEX: inside YYt_update_local: R_y.*R_y_mask');

Pwrt_Ylk=-2*(transpose(R_y)*Y_lk)-2*(R_y*Y_lk)+4*(Y_lk*transpose(Y_lk)*Y_lk);              

Y_lk=Y_lk*exp(-aeta*Pwrt_Ylk);

end


%find the indices of the L_target=3 most dominant diagonal elements.
%compute their average amplitude/energy
%Set R_y at these indices equal to the average: a crude normalization. 
function [R_y]=YYt_update_sort_diagonal( R_y,  L_static)

L_target=3;

diagvec=diag(R_y);

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

%use the 2nd column of the info column matrix to index and set the
%avg_energy value. 
diagvec(squeeze(Ltargetx2_info_colmatrix(:, 2)))=avg_energy;

%Set R_y at these indices equal to the average: a crude normalization. 
R_y(logical(eye(L_static))) = diagvec;

end