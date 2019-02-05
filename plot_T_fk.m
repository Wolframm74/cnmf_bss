function plot_T_fk(T_fk, close_all_flag)

if (close_all_flag)
    %close all;

    close_figures_local();
    
end
        figure,
        %subplot(3,1,1)
        imagesc(T_fk);

        title('MEX: T_fk colormap');

%         figure, hist(T_fk)
        
% T_fk        
%         
%  waitforbuttonpress;         
        
end

function close_figures_local()

h=findobj('type','figure');
n=length(h);

if (n>30)

    close all; 
    
end

end