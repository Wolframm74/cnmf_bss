function [Y_lk]=orthogonalize_Y_lk_square(Y_lk)

display('INSIDE ORTHOGONALIZE Y_LK SQUARE')

[L_static, K_static]=size(Y_lk);

% if (K_static==L_static)
    
    Y_lk=orthogonalize_Y_lk_squaremat_soft(Y_lk);
    
% end

end
%immediately after returning mex application should 

function [Y_lk]= orthogonalize_Y_lk_squaremat(Y_lk)

[U_y, S_y, V_y]=svd(Y_lk);

Y_lk=U_y*transpose(V_y);

%search for indices w/ negative elements
X=Y_lk<=0;

%zero these indices.
Y_lk(X)=zeros();

end

function [Y_lk]= orthogonalize_Y_lk_squaremat_soft(Y_lk)

[U_y, S_y, V_y]=svd(Y_lk);

%Orthogonalize just the U_y matrix!
U_y=orthogonalize_Y_lk_squaremat(U_y);

%Keep the signular values as part of the SVD representation!
Y_lk=U_y*S_y*transpose(V_y);

%search for indices w/ negative elements
X=Y_lk<=0;

%zero these indices.
Y_lk(X)=zeros();

end
