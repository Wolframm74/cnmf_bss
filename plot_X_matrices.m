function plot_X_matrices(Xhat_outtensor_cx_FLM, Xhat_outtensor_real_FLM, Xhat_outtensor_cx_FKM, Xhat_outtensor_real_FKM)

        figure,
        subplot(2,1,1)
        imagesc(abs(Xhat_outtensor_cx_FLM(:,:,1)));              
        title('MLref: Xhat_outtensor_cx_FLM(:,:,1) colormap');
        
        subplot(2,1,2)
        imagesc(abs(Xhat_outtensor_cx_FLM(:,:,2)));              
        title('MLref: Xhat_outtensor_cx_FLM(:,:,2) colormap');        
        
        figure,
        subplot(2,1,1)
        imagesc(Xhat_outtensor_real_FLM(:,:,1));              
        title('MLref: Xhat_outtensor_real_FLM(:,:,1) colormap');
        
        subplot(2,1,2)
        imagesc(Xhat_outtensor_real_FLM(:,:,2));              
        title('MLref: Xhat_outtensor_real_FLM(:,:,2) colormap');        
        
        figure,
        subplot(2,1,1)
        imagesc(abs(Xhat_outtensor_cx_FKM(:,:,1)));              
        title('MLref: Xhat_outtensor_cx_FKM(:,:,1) colormap');
        
        subplot(2,1,2)
        imagesc(abs(Xhat_outtensor_cx_FKM(:,:,2)));              
        title('MLref: Xhat_outtensor_cx_FKM(:,:,2) colormap');        
        
        figure,
        subplot(2,1,1)
        imagesc(Xhat_outtensor_real_FKM(:,:,1));              
        title('MLref: Xhat_outtensor_real_FKM(:,:,1) colormap');
        
        subplot(2,1,2)
        imagesc(Xhat_outtensor_real_FKM(:,:,2));              
        title('MLref: Xhat_outtensor_real_FKM(:,:,2) colormap');                
end