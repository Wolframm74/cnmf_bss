function [Z_ol]=orthogonalize_Z_ol_skinny(Z_ol)

L_target=3;

display('INSIDE ORTHOGONALIZE Z_ol SQUARE')

[O_static, L_static]=size(Z_ol);
   
[U_z, S_z, V_z]=svd(Z_ol);

princip_sing_values_sum=0;

for i=1:L_static
    
    if i<=L_target
        
        princip_sing_values_sum=princip_sing_values_sum+abs(S_z(i,i));
        
    end
    
    %Try zeroing singular values greater than the desired number of
    %sources. We would like Z_ol to have column rank=L_target; 
    if i>=L_target
        
        %S_z(i,i)=S_z(i,i)*exp(-((i-L_target))/1.5);
        
        S_z(i,i)=S_z(i,i)*((-1/(L_static-L_target))*(i-L_target)+1);
        
        %instead of zeroing, scale the "residual" singular values.
        
        %test code:
        %x=0:24
        %figure, plot(x, exp(-x/1.5))
        
        
    end
        
end

for i=1:L_target
    
    S_z(i,i)=sign(S_z(i,i))*(princip_sing_values_sum/L_target);
    
end

%Reconstruct Z_ol. 
Z_ol=U_z*S_z*transpose(V_z);

    %since we expect that L_static<O_static (for a skinny matrix, pass in
    %L_static, since the number of singular values should equal at most th
    %e less between L_static and O_static. 
    Z_ol=orthogonalize_Z_ol_soft_wrapper(Z_ol, L_static);

end
%immediately after returning mex application should 

%don't zero the elements of U_z. can allow it to be negative
%V_z should be of dimensions LxL. 
function [V_z]= orthogonalize_squaremat_hard(V_z)

[U, S, V]=svd(V_z);

%this will forcefully impose that the output matrix V_z, should you call svd() again
%on it now: I think will have S as the identity matrix!
V_z=U*transpose(V);

%search for indices w/ negative elements
%X=Z_ol<=0;

%zero these indices.
%Z_ol(X)=zeros();
%don't zero negative elements of the orthogonal output matrix

end

function [Z_ol]= orthogonalize_Z_ol_soft_wrapper(Z_ol, L_static)

[U_z, S_z, V_z]=svd(Z_ol);

%Orthogonalize just the U_z matrix!
V_z=orthogonalize_squaremat_hard(V_z);

%Keep the signular values as part of the SVD representation!
Z_ol=U_z*S_z*transpose(V_z);

%Though the V_z matrix is allowed to have negative values: the output
%matrix Z_ol should not be allowed to!
%search for indices w/ negative elements
X=Z_ol<=0;

%zero these indices.
Z_ol(X)=zeros();

end
