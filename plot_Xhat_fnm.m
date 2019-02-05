function plot_Xhat_fnm(Xhat_fnm)

        figure,
        subplot(2,1,1)
        imagesc(abs(Xhat_fnm(:,:,1)));
        title('MEX: Xhat_fnm colormap');
                
        subplot(2,1,2)
        imagesc(abs(Xhat_fnm(:,:,2)));

end