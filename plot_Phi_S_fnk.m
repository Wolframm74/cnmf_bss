%total amount of figures will be figure_num*subplot_num
function plot_Phi_S_fnk(expj_Phi_S_fkn, figure_num, subplot_num)

figure_num=4;

subplot_num=6;

%K=30=5*6

outer_iter=0;

inner_iter=0;

plot_Phi_S_fnk_2(expj_Phi_S_fkn, 1, 4);

%figure,

% for outer_iter=1:figure_num
% 
%     figure, 
%     
%     for inner_iter=1:subplot_num
% 
%        index=outer_iter*inner_iter;
%         
%        subplot(5,6, index) 
%        subplot(1,subplot_num, inner_iter) 
%        
%        imagesc(angle(squeeze(expj_Phi_S_fkn(:,index,:))));       
%              
%     end
%         
% end
% 
% title('MEX: Phi_S_fnk(:,:,k=1:30) colormap');  

%waitforbuttonpress         

%         figure, 
%         %subplot(2,1,1)
%         imagesc(abs(Phi_S_fnk(:,:,1)));
%         
%         title('MEX: Phi_S_fnk(:,:,1) colormap');                     
        
end

function plot_Phi_S_fnk_2(expj_Phi_S_fkn, k_start, k_end)

for k_index=k_start:k_end
    
    figure,
    imagesc(angle(squeeze(expj_Phi_S_fkn(:,k_index,:))));      
    title([ num2str(k_index), 'th component channel, phase']);
    
end


end