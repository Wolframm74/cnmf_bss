function [V_nk]=VVt_update(V_nk)

[N_static, K_static]=size(V_nk);

V_kn=transpose(V_nk);

R_v=V_kn*transpose(V_kn);

R_v_indices=find(R_v>25);

R_v(R_v_indices)=25;

R_v_mask=zeros(K_static, K_static);

R_v_mask=diag(ones(K_static, 1), 0);

% R_v_mask=R_v_mask+diag(0.5*ones(K_static-1, 1), 1);
% 
% R_v_mask=R_v_mask+diag(0.5*ones(K_static-1, 1), -1);
% 
% R_v_mask=R_v_mask+diag(0.25*ones(K_static-2, 1), 2);
% 
% R_v_mask=R_v_mask+diag(0.25*ones(K_static-2, 1), -2);
% 
% R_v_mask=R_v_mask+tril(0.1*ones(K_static,K_static), -3);
% 
% R_v_mask=R_v_mask+triu(0.1*ones(K_static,K_static), 3);

for i=1:(K_static-2)
    
% R_v_mask=R_v_mask+diag((1-i*(1/K_static))*ones(K_static-i, 1), i);
% 
% R_v_mask=R_v_mask+diag((1-i*(1/K_static))*ones(K_static-i, 1), -i);

R_v_mask=R_v_mask+diag(exp(-i/3)*ones(K_static-i, 1), i);

R_v_mask=R_v_mask+diag(exp(-i/3)*ones(K_static-i, 1), -i);

%{
    diagvec=exp(-i/(8))*ones(K_static-i, 1);
    %diagvec=(1-i*(1/K_static))*ones(K_static-i, 1);
    
    for k=1:(K_static-i)
        
        diagvec(k)=((k-(K_static-i+1)/2)^2)*diagvec(k);
        
    end

start_index=((1)-(K_static-i+1)/2);
offset=1-2*(i-1)*(1/(K_static-2));    

if ((offset>1)||(offset<0))
    
    offset=0;
    
end

alpha=(1-offset)/(start_index^2);

%normalize the diagvec by its first element's amplitude
diagvec=diagvec*alpha+offset*ones(numel(diagvec), 1);
    
R_v_mask=R_v_mask+diag(diagvec, i);

R_v_mask=R_v_mask+diag(diagvec, -i);

%}
end

        figure, 
        imagesc(V_kn);
        title('MEX: V_kn: before');
        
        figure,
        %subplot(3,1,1)
        imagesc(R_v);

        title('MEX: VVt_update: before');
        
        figure,
        %subplot(3,1,1)
        imagesc(R_v_mask);

        title('MEX: R_v_mask');

V_kn=VVt_update_local(V_kn, R_v.*R_v_mask, K_static, N_static);        
R_v=V_kn*transpose(V_kn);

        figure, 
        imagesc(V_kn);
        title('MEX: V_kn: after');

        figure,
        %subplot(3,1,1)
        imagesc(R_v);

        title('MEX: VVt_update: after');
        
        V_nk=transpose(V_kn);

        
end

function [V_kn]=VVt_update_local(V_kn, R_v, K_static, N_static)

gamma=10;
%gamma=50;

figure, 
imagesc(R_v);
title('MEX: VVt_update: R_v before normalization');

%[R_v]=VVt_update_normalize_diagonal( R_v,  K_static);

figure, 
imagesc(R_v);
title('MEX: VVt_update: R_v after normalization');

        figure,
        %subplot(3,1,1)
        imagesc(R_v);
        title('MEX: inside VVt_update_local: R_v.*R_v_mask');

onesmat_KN=ones(K_static,N_static);

V_kn=V_kn.*((1-(1/gamma))*onesmat_KN+(1/(2*gamma))*(transpose(R_v)*V_kn+R_v*V_kn)./(V_kn*transpose(V_kn)*V_kn));

end

function [R_v]=VVt_update_normalize_diagonal( R_v,  K_static)

aeta=0.1;

diagvec=diag(R_v);

figure, 
plot(diagvec);
title('diagvec before');

mean_diagvec=(1/K_static)*sum(diagvec);

Partialvec=2*(diagvec-mean_diagvec*ones(K_static, 1));

diagvec=diagvec.*exp(-aeta*Partialvec);

figure, 
plot(diagvec);
title('diagvec after gradient update');

%compute the mean again
mean_diagvec=(1/K_static)*sum(diagvec);

[diagvec]=elementwise_scaling_diagvec(diagvec, mean_diagvec, K_static);

figure, 
plot(diagvec);
title('diagvec after');

%Set R_y at these indices equal to the average: a crude normalization. 
R_v(logical(eye(K_static))) = diagvec;

end

function [diagvec]=elementwise_scaling_diagvec(diagvec, mean_diagvec, K_static)

mu=(1/12)*mean_diagvec;

sigma=mu/4;

beta=8;

elementwise_scaling_vector=0.5*(1+erf((diagvec-mu)/(sigma*sqrt(2))));

figure, 
plot(elementwise_scaling_vector);
title('elementwise_scaling_vector');

diagvec=elementwise_scaling_vector.*diagvec;

%give a boost to the middle area to prevent small to medium clusters from
%dissapearing.

alpha=mean_diagvec/2;

elementwise_scaling_vector=0.3*exp(-((abs(diagvec-mean_diagvec)/alpha).^beta));
%elementwise_scaling_vector=0.2*exp(-((diagvec-mean_diagvec)/2*sigma*sigma).^2);

elementwise_scaling_vector=elementwise_scaling_vector+0.9*ones(K_static, 1);

diagvec=elementwise_scaling_vector.*diagvec;

end