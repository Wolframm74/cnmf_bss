function [T_fk]=TtT_update(T_fk)

display('inside TtT_update');

[F_static, K_static]=size(T_fk);

R_t=transpose(T_fk)*T_fk;

R_t_indices=find(R_t>10000000);

R_t(R_t_indices)=10000000;

R_t_mask=zeros(K_static, K_static);

R_t_mask=diag(ones(K_static, 1), 0);

for i=1:(K_static-2)

R_t_mask=R_t_mask+diag(exp(-i/2)*ones(K_static-i, 1), i);

R_t_mask=R_t_mask+diag(exp(-i/2)*ones(K_static-i, 1), -i);

end

        figure, 
        imagesc(T_fk);
        title('MEX: T_fk: before');

        figure,
        imagesc(R_t);        
        title('MEX: TtT_update: before');
        
        figure,
        imagesc(R_t_mask);
        title('MEX: R_t_mask');        
        
T_fk=TtT_update_local(T_fk, R_t.*R_t_mask, F_static, K_static);        

R_t=transpose(T_fk)*T_fk;

        figure, 
        imagesc(T_fk);
        title('MEX: T_fk: after');

        figure,
        imagesc(R_t);        
        title('MEX: TtT_update: after');
        

end

function [T_fk]=TtT_update_local(T_fk, R_t, F_static, K_static)

gamma=50;

figure, 
imagesc(R_t);
title('MEX: inside TtT_update_local: R_t.*R_t_mask');

onesmat_FK=ones(F_static, K_static);

T_fk=T_fk.*((1-(1/gamma))*onesmat_FK+(1/(2*gamma))*(T_fk*R_t+T_fk*transpose(R_t))./(T_fk*transpose(T_fk)*T_fk));


end