function [L_value_2nd_costfun]=plot_Xhat_Xtilde_fnm(Xhat_fnm, Xtilde_fnm)

L_value_2nd_costfun=zeros(3,1);

dummy_value=0;

[F_static, N_static, M_static]=size(Xhat_fnm);

Xtilde_fnm_phase_normalized=zeros(F_static, N_static);

        figure,
        subplot(2,1,1)
        imagesc(abs(Xtilde_fnm(:,:,1)));            
        title('MEX: Xtilde_fnm colormap');
                
        subplot(2,1,2)
        imagesc(abs(Xtilde_fnm(:,:,2)));

        figure,
        subplot(2,1,1)
        imagesc(abs(Xhat_fnm(:,:,1)));
        title('MEX: Xhat_fnm colormap');
                
        subplot(2,1,2)
        imagesc(abs(Xhat_fnm(:,:,2)));

        figure, 
        imagesc(angle(Xtilde_fnm(:,:,1)));
        title('MEX: imagesc(angle(Xtilde_fnm(:,:,1)))');        
        
        figure, 
        imagesc(angle(Xtilde_fnm(:,:,2)));
        title('MEX: imagesc(angle(Xtilde_fnm(:,:,2)))');
        
        figure, 
        imagesc(angle(Xtilde_fnm(:,:,1))-angle(Xtilde_fnm(:,:,2)));
        title('MEX: imagesc(angle(Xtilde_fnm(:,:,1))-angle(Xtilde_fnm(:,:,2)))');
        
        figure, 
        imagesc(angle(Xhat_fnm(:,:,1)));
        title('MEX: imagesc(angle(Xhat_fnm(:,:,1)))');        
        
        figure, 
        imagesc(angle(Xhat_fnm(:,:,2)));
        title('MEX: imagesc(angle(Xhat_fnm(:,:,2)))');        
        
        figure, 
        imagesc(angle(Xhat_fnm(:,:,1))-angle(Xhat_fnm(:,:,2)));    
        title('MEX: imagesc(angle(Xhat_fnm(:,:,1))-angle(Xhat_fnm(:,:,2)))');
        
%         [Xhat_fnm_phase_normalized, expj_Xtilde_fnm]=frobenius_norm_X_fnm_phase_normalization(Xhat_fnm, F_static, N_static);
%         
%         figure, 
%         imagesc(angle(Xhat_fnm_phase_normalized(:,:,2)));
%         title('MEX: imagesc(angle(Xhat_fnm_phase_normalized(:,:,2)))');    

        Xhat_fnm_phase_normalized=squeeze(Xhat_fnm(:, :, 2)).*conj(squeeze(Xhat_fnm(:, :, 1)));

        figure, 
        imagesc(angle(Xhat_fnm_phase_normalized));
        %title('MEX: imagesc(angle(Xhat_fnm_phase_normalized))');            
        title('Interchannel Phase Difference Quantity: arg( X_b(f,n) conj (X_a(f,n) ) ');            
        ylabel('Frequency');
        xlabel('Time Activation');

        Xtilde_fnm_phase_normalized=squeeze(Xtilde_fnm(:, :, 2)).*conj(squeeze(Xtilde_fnm(:, :, 1)));

        figure, 
        imagesc(angle(Xtilde_fnm_phase_normalized));
        %title('MEX: imagesc(angle(Xtilde_fnm_phase_normalized))');            
        title('Interchannel Phase Difference Quantity: arg( X_b(f,n) conj (X_a(f,n) ) ');    
        ylabel('Frequency');
        xlabel('Time Activation');
        
        [dummy_value]=compute_L_2nd_local(Xhat_fnm, Xtilde_fnm, 2, 1);
        
        L_value_2nd_costfun(1)=dummy_value;
        
        display(L_value_2nd_costfun)
        
end

function [Xtilde_fnm_phase_normalized, expj_Xtilde_fnm]=frobenius_norm_X_fnm_phase_normalization(Xtilde_fnm, F_static, N_static)

%Strip Xtilde_fnm of its absolute value and save its phase at each FxN bin.
expj_Xtilde_fnm=Xtilde_fnm./abs(Xtilde_fnm);

channel_m_equals_1_index=1; 

for f_index=1:F_static
    
    for n_index=1:N_static

        %frob_norm_mat_FN(f_index, n_index)=sqrt(sum(abs(squeeze(Xtilde_fnm(f_index, n_index, :)).*abs(squeeze(Xtilde_fnm(f_index, n_index, :))))));        
        
        %first normalize by its frobenius norm
        %Xtilde_fnm_normalized(f_index, n_index, :)=squeeze(Xtilde_fnm(f_index, n_index, :))/frob_norm_mat_FN(f_index, n_index);
        
        %remove the phase of the first channel by multiplying by the
        %conjugate
        Xtilde_fnm_phase_normalized(f_index, n_index, :)=Xtilde_fnm(f_index, n_index, :)*conj(expj_Xtilde_fnm(f_index, n_index, 1));
        
        %can restore the phase of the first channel by multiplying
        %Xtilde_fnm_nomalized by expj_Xtilde_fnm(f_index, n_index, 1), when
        %you return from NMF. 
        
    end
    
end

end 

function [L_value_2nd_costfun]=compute_L_2nd_local(Xhat_fnm, Xtilde_fnm, b_value, a_value)

L_value_2nd_costfun=sum(sum((squeeze(Xtilde_fnm(:, :, b_value)).*(conj(squeeze(Xtilde_fnm(:, :, a_value))))-squeeze(Xhat_fnm(:, :, b_value)).*(conj(squeeze(Xhat_fnm(:, :, a_value)))) ).^2)); 

end

