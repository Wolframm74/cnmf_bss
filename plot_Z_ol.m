function [Z_ol]=plot_Z_ol(Z_ol)

        figure,
        %subplot(3,1,1)
        imagesc(Z_ol);        
        
        title('MEX: Z_ol colormap');

        figure,
        %subplot(3,1,1)
        waterfall(transpose(Z_ol));        
        
        title('MEX: waterfall(Z_ol) ');
        
%         figure,
%         hist(Z_ol)

% Z_ol      
%         
% waitforbuttonpress        
        
end