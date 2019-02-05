function [Z_ol]=ZtZ_update(Z_ol)

[O_static, L_static]=size(Z_ol);

[U_z S_z V_z]=svd(Z_ol);

for l_iter=1:L_static
    
S_z(l_iter, l_iter)=1;   
    
end

end