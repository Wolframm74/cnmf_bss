function [V_nk]=orthogonalize_V_nk(V_nk)

V_kn=transpose(V_nk);

[N_static, K_static]=size(V_nk);

display('INSIDE orthogonalize_V_nk')
    
    V_kn=orthogonalize_Y_lk_square(V_kn);
    
V_nk=transpose(V_kn);

end
%immediately after returning mex application should 

% function [V_nk]= orthogonalize_Vkn_fatmat(V_kn, K_static)
% 
% [U_y, S_y, V_y]=svd(V_kn);
% 
% for i=1:K_static
%     
%     if (S_y(i,i)~=0)
%         
%         S_y(i,i)=1;
%         
%     end
%     
% end
% 
% V_kn=U_y*S_y*transpose(V_y);
% 
% search for indices w/ negative elements
% X=V_kn<=0;
% 
% zero these indices
% V_kn(X)=zeros();
% 
% V_nk=transpose(V_kn);
% 
% end
