function [Z_ol]=ZtZ_update(Z_ol)

R_z=transpose(Z_ol)*Z_ol;

        figure,
        %subplot(3,1,1)
        imagesc(R_z);        
        
        title('MEX: ZtZ_update');
        
end