%potential output arguments: outtensor_2nd_real_FKN

function Phi_S_unit_test(expj_Phi_S_fkn, T_fk, V_nk, Xhat_low_fnm, E_conj_fnm, Xhat_outtensor_real_MKF, Xhat_outtensor_cx_MKF, Y_lk, Z_ol, W_fom, expj_Phi_W_fom, dim_mat)

%save('Phi_S_unit_test.mat');

%waitforbuttonpress         

M_static=dim_mat(1);
F_static=dim_mat(2);
N_static=dim_mat(3);
K_static=dim_mat(4);
L_static=dim_mat(5);
O_static=dim_mat(6); 

outtensor_2nd_real_FKN=zeros(F_static, K_static, N_static);
Xhat_low_fmn_local=zeros(F_static, M_static, N_static);

for n_index=1:N_static

    for m_index=1:M_static

       Xhat_low_fmn_local(:, m_index, n_index)=squeeze(Xhat_low_fnm(:, n_index, m_index));
       
       E_conj_fmn_local(:, m_index, n_index)=squeeze(E_conj_fnm(:, n_index, m_index));

    end
    
    for f_index=1:F_static
        
    outtensor_2nd_real_FKN(f_index, :, n_index)=squeeze(Xhat_low_fmn_local(f_index,:, n_index))*squeeze(Xhat_outtensor_real_MKF(:,:,f_index));

    outtensor_2nd_cx_FKN(f_index, :, n_index)=squeeze(E_conj_fmn_local(f_index,:, n_index))*squeeze(Xhat_outtensor_cx_MKF(:,:,f_index));

    tensor_dummy_cx_FKN(f_index, :, n_index)=squeeze(outtensor_2nd_cx_FKN(f_index, :, n_index)).*squeeze(T_fk(f_index,:)).*squeeze(V_nk(n_index,:));

    outtensor_target_cx_FKN(f_index, :, n_index)=conj(squeeze(expj_Phi_S_fkn(f_index, :, n_index))).*squeeze(T_fk(f_index,:)).*squeeze(V_nk(n_index,:));        

    outtensor_target_cx_FKN(f_index, :, n_index)=squeeze(outtensor_target_cx_FKN(f_index, :, n_index)).*squeeze(outtensor_2nd_real_FKN(f_index, :, n_index));
    
    outtensor_target_cx_FKN(f_index, :, n_index)=squeeze(outtensor_target_cx_FKN(f_index, :, n_index))+squeeze(tensor_dummy_cx_FKN(f_index, :, n_index));
    
    end

end

for f=1:F_static
   
    for n=1:N_static
       
        for k=1:K_static
        
            if (abs(outtensor_target_cx_FKN(f,k,n))==0)
               
                outtensor_target_cx_FKN(f,k,n)=-1i;
                
            end
            
        end
                
    end
    
end

phase_shift_tensor_cx_FKN=(-1)*ones(F_static, K_static, N_static);

reference_phase_tensor_cx_FKN=(1i)*ones(F_static, K_static, N_static);

absval_tensor_FKN=abs(outtensor_target_cx_FKN);

expj_Phi_S_fkn=(reference_phase_tensor_cx_FKN./(outtensor_target_cx_FKN./absval_tensor_FKN)).*phase_shift_tensor_cx_FKN;

%save('Phi_S_unit_test.mat');

display('Got to checkpoint: Phi_S_unit_test')

plot_Phi_S_fnk_local(expj_Phi_S_fkn);

end

function plot_Phi_S_fnk_local(expj_Phi_S_fkn)

%K=30=5*6

outer_iter=0;

inner_iter=0;

%figure,

for outer_iter=1:5

    figure, 
    
    for inner_iter=1:6

       index=outer_iter*inner_iter;
        
       %subplot(5,6, index) 
       subplot(1,6, inner_iter) 
       
       imagesc(angle(squeeze(expj_Phi_S_fkn(:,index,:))));       
             
    end
        
end

title('MEX: Phi_S_fnk(:,:,k=1:30) colormap');  

waitforbuttonpress         

%         figure, 
%         %subplot(2,1,1)
%         imagesc(abs(Phi_S_fnk(:,:,1)));
%         
%         title('MEX: Phi_S_fnk(:,:,1) colormap');                     
        
end