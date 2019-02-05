function plot_W_fom(W_fom)

        figure,
        subplot(2,1,1)
        imagesc(abs(W_fom(:,:,1)));
        grid on
        title('MEX: W_fom colormap');
                        
        subplot(2,1,2)
        imagesc(abs(W_fom(:,:,2)));
        
%         figure, hist(W_fom)
        
% W_fom        
%         
% waitforbuttonpress                 
        
end