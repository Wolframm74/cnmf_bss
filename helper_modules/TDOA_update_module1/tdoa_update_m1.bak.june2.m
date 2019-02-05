function [Z_ol, W_fom] = tdoa_update_m1(Xtilde_fnm, Z_ol, Y_lk, T_fk, V_nk, W_fom, expj_Phi_W_fom, expj_Phi_S_fkn)

plotfun_local_Z_W(Z_ol, W_fom, 1);

[N_static, K_static]=size(V_nk);

[F_static, K_static]=size(T_fk);

[L_static, K_static]=size(Y_lk);

[F_static, O_static, M_static]=size(W_fom);

Y_est_FxNxL=zeros(F_static, N_static, L_static);

tempmat_FN=zeros(F_static, N_static);

for l_index=1:L_static
    
    tempmat_FN=zeros(F_static, N_static);
    
    for k_index=1:K_static
        
        tempmat_FN=tempmat_FN+Y_lk(l_index, k_index)*kron( T_fk(:, k_index ) , transpose( V_nk(:, k_index) ) ).*squeeze(expj_Phi_S_fkn(:, k_index, :));        
        
    end
    
    Y_est_FxNxL(:, :, l_index)=tempmat_FN;
    
end

%% Compute A

A_MLF=zeros(M_static, L_static, F_static);

MxL_dummy_mat1=zeros(M_static, L_static);
LxL_dummy_mat1=zeros(L_static, L_static);

for f_index=1:F_static
    
    MxL_dummy_mat1=zeros(M_static, L_static);
    LxL_dummy_mat1=zeros(L_static, L_static);
    
    for n_index=1:N_static
        
        MxL_dummy_mat1=MxL_dummy_mat1+squeeze(Xtilde_fnm(f_index, n_index, :))*(squeeze(Y_est_FxNxL(f_index, n_index, :))');
        
        LxL_dummy_mat1=LxL_dummy_mat1+squeeze(Y_est_FxNxL(f_index, n_index, :))*(squeeze(Y_est_FxNxL(f_index, n_index, :))');       
        
    end
     
    A_MLF(:, :, f_index)=((1/N_static)*MxL_dummy_mat1)*pinv((1/N_static)*LxL_dummy_mat1);   
    
end

wlen=4*512;
h=wlen/4;
%h=wlen;
nfft=wlen;
Fs=16000;

J_value=1;
N_stft=nfft; 

rhat_LMF=zeros(L_static, M_static, F_static);
b_LMF=zeros(L_static, M_static, F_static);

freq=zeros(F_static, 1);

for f_index=2:F_static
    
    freq(f_index)=((f_index-1)*Fs)/(N_stft);
    
    for l_index=1:L_static
        
        for m_index=1:M_static    
            
        rhat_LMF(l_index, m_index, f_index)=-angle( A_MLF(m_index, l_index, f_index) / A_MLF(J_value, l_index, f_index ) )/(2*pi*freq(f_index));    
        
        b_LMF(l_index, m_index, f_index)=exp(-1i*2*pi*freq(f_index)*rhat_LMF(l_index, m_index, f_index));
        
        end
    
    end
        
end




b_LMF_target=zeros(L_static, M_static, F_static);

bhat_ML=zeros(M_static, L_static);

dummy_Mx1=zeros(M_static, 1);

max_TDOA=0.0015;

for l_index=1:L_static
    
    dummy_Mx1=zeros(M_static, 1);
    
    for f_index=1:F_static
        
        dummy_Mx1=dummy_Mx1+transpose(squeeze(b_LMF(l_index, :, f_index)));
        
    end
    
    bhat_ML(:, l_index)=(1/F_static)*dummy_Mx1;
    
    %for m=2 clip the max TDOA value negative or positive.
    
    if (abs(bhat_ML(2, l_index))>max_TDOA)
        
        bhat_ML(2, l_index)=sign(bhat_ML(2, l_index))*max_TDOA;
        
    end
    
    b_LMF_target(l_index, :, :)=kron( bhat_ML(:, l_index), ones(1, F_static) );
    
end

%% compute a few averages to detect if things have converged down enough
avg1=sum(abs(bhat_ML(:,1)));

avg2=sum(abs(bhat_ML(:,2)));

avg3=sum(abs(bhat_ML(:,3)));

avg_total=avg1+avg2+avg3;

%if (avg_total<0.008*3) %only do the update if we have converged down enough

%% Scatter plots for f=1:F

figure,
scatter(squeeze(rhat_LMF(1, 1, :)), squeeze(rhat_LMF(1, 2, :)));

figure,
scatter(squeeze(rhat_LMF(2, 1, :)), squeeze(rhat_LMF(2, 2, :)));

figure,
scatter(squeeze(rhat_LMF(3, 1, :)), squeeze(rhat_LMF(3, 2, :)));

rhat_FLM=zeros(F_static, L_static, M_static);

for f_index=1:F_static
    
    rhat_FLM(f_index, :, :)=squeeze(rhat_LMF(:, :, f_index));
    
end

rhat_tallmat=reshape(rhat_FLM, F_static*L_static, M_static);

figure, 
scatter(rhat_tallmat(:, 1), rhat_tallmat(:,2));
    
%% Scatter plots f=1:F averaged out.

figure, 
scatter( bhat_ML(1,1), bhat_ML(2,1) );

figure, 
scatter( bhat_ML(1,2), bhat_ML(2,2) );

figure, 
scatter( bhat_ML(1,3), bhat_ML(2,3) );

%% compute h_flm

h_flm=zeros(F_static, L_static, M_static);

dummy_Mx1=zeros(M_static, 1);
    
for f_index=1:F_static

    for l_index=1:L_static
        
        dummy_Mx1=zeros(M_static, 1);
        
        for o_index=1:O_static
            
        dummy_Mx1=dummy_Mx1+Z_ol(o_index, l_index)*squeeze(W_fom(f_index, o_index, :));
            
        end
        
        h_flm(f_index, l_index, :)=dummy_Mx1;
    
    end
    
end

%% call the update functions

[W_fom] = W_fom_update_local(b_LMF_target, h_flm, Z_ol, W_fom, expj_Phi_W_fom, F_static, O_static, M_static, L_static);

[Z_ol] = Z_ol_update_local(b_LMF_target, h_flm, Z_ol, W_fom, expj_Phi_W_fom, O_static, L_static, F_static);

%% plot check
plotfun_local_Z_W(Z_ol, W_fom, 2);

%end %end if

end

function [W_fom] = W_fom_update_local(b_LMF_target, h_flm, Z_ol, W_fom, expj_Phi_W_fom, F_static, O_static, M_static, L_static)

partial_d_wrt_Wfom_1=zeros(F_static, O_static);

partial_d_wrt_Wfom_2=zeros(F_static, O_static);

J_value=1;

dummy_Lx1_1=zeros(L_static, 1);
dummy_Lx1_2=zeros(L_static, 1);

%for m_index=2:M_static
m_index=2;
    
    for o_index=1:O_static
     
        for f_index=1:F_static

            %dummy_Lx1_1=zeros(L_static, 1);
            
            dummy_Lx1_1=conj(squeeze(b_LMF_target(:, m_index, f_index)));
            
            dummy_Lx1_1=dummy_Lx1_1./transpose(squeeze(h_flm(f_index, :, m_index)));
            
            dummy_Lx1_1=expj_Phi_W_fom(f_index, o_index, m_index)*dummy_Lx1_1;
            
            dummy_Lx1_1=dummy_Lx1_1.*transpose(Z_ol(o_index, :));
            
            %partial_d_wrt_Wfom_1(f_index, o_index, m_index)=sum(-2*real(dummy_Lx1_1));
            partial_d_wrt_Wfom_1(f_index, o_index)=sum(-2*real(dummy_Lx1_1));
            
            dummy_Lx1_2=transpose(squeeze(h_flm(f_index, :, m_index)));
            
            dummy_Lx1_2=dummy_Lx1_2./transpose(squeeze(h_flm(f_index, :, J_value)));
            
            dummy_Lx1_2=2*dummy_Lx1_2./abs(dummy_Lx1_2);
            
            dummy_Lx1_2=dummy_Lx1_2./transpose(squeeze(h_flm(f_index, :, J_value)));
            
            dummy_Lx1_2=expj_Phi_W_fom(f_index, o_index, m_index)*dummy_Lx1_2;
            
            dummy_Lx1_2=dummy_Lx1_2.*transpose(Z_ol(o_index, :));
            
            dummy_Lx1_2=real(dummy_Lx1_2);
            
            %partial_d_wrt_Wfom_2(f_index, o_index, m_index)=sum(dummy_Lx1_2);
            partial_d_wrt_Wfom_2(f_index, o_index)=sum(real(dummy_Lx1_2));
            
        end
        
    end
    
%end
    
%do gradient descent here on W_fom
%W_fom=> but only update the m=2th slice.

%gamma=50;
gamma=10;

epsilon_mat=zeros(F_static, O_static);

%elementwise division by the positive part of the partial derivative.
epsilon_mat=(1/gamma)*W_fom(:, :, m_index)./partial_d_wrt_Wfom_2;

aeta=0.2;
%aeta=1;

%W_fom(:, :, m_index)=squeeze(W_fom(:, :, m_index)).*exp(-aeta*(partial_d_wrt_Wfom_1+partial_d_wrt_Wfom_2));

W_fom(:, :, m_index)=W_fom(:, :, m_index)-epsilon_mat.*(partial_d_wrt_Wfom_1+partial_d_wrt_Wfom_2);

end

function [Z_ol] = Z_ol_update_local(b_LMF_target, h_flm, Z_ol, W_fom, expj_Phi_W_fom, O_static, L_static, F_static)

partial_d_wrt_Zol_1=zeros(O_static, L_static);

partial_d_wrt_Zol_2=zeros(O_static, L_static);

m_index=2;

J_value=1;

dummy_Fx1_1=zeros(F_static, 1);
dummy_Fx1_2=zeros(F_static, 1);
dummy_Fx1_3=zeros(F_static, 1);

for o_index=1:O_static
     
        for l_index=1:L_static
            
            dummy_Fx1_1=conj(squeeze(b_LMF_target(l_index, m_index, :)));
            
            dummy_Fx1_3=( squeeze(h_flm(:,l_index, J_value)).*squeeze(W_fom(:, o_index, m_index)).*squeeze(expj_Phi_W_fom(:, o_index, m_index)) - squeeze(h_flm(:, l_index, m_index)).*squeeze(W_fom(:, o_index, J_value)).*squeeze(expj_Phi_W_fom(:, o_index, J_value)) );
            
            dummy_Fx1_3=dummy_Fx1_3./( squeeze(h_flm(:, l_index, J_value)).^2 );
            
            dummy_Fx1_1=dummy_Fx1_1.*dummy_Fx1_3;
            
            partial_d_wrt_Zol_1(o_index, l_index)=sum(-2*real(dummy_Fx1_1));
            
            dummy_Fx1_2=squeeze(h_flm(:, l_index, m_index))./squeeze(h_flm(:, l_index, J_value));
            
            dummy_Fx1_2=2*dummy_Fx1_2./abs(dummy_Fx1_2);
            
            dummy_Fx1_2=dummy_Fx1_2.*dummy_Fx1_3;
            
            partial_d_wrt_Zol_2(o_index, l_index)=sum(real(dummy_Fx1_2));
            
        end
        
end

%do gradient descent

%gamma=50;
gamma=10;

epsilon_mat=zeros(O_static, L_static);

%elementwise division by the positive part of the partial derivative.
epsilon_mat=(1/gamma)*Z_ol./partial_d_wrt_Zol_2;

figure, imagesc(partial_d_wrt_Zol_1);

figure, imagesc(partial_d_wrt_Zol_2);

figure, imagesc(partial_d_wrt_Zol_1 + partial_d_wrt_Zol_2);

%Z_ol=Z_ol.*exp(-aeta*(partial_d_wrt_Zol_1+partial_d_wrt_Zol_2));
Z_ol=Z_ol-epsilon_mat.*(partial_d_wrt_Zol_1+partial_d_wrt_Zol_2);

end

function plotfun_local_Z_W(Z_ol, W_fom, input_flag)

if (input_flag<1.5)
    
    string1=' tdoa_update_m1: before ';
    
end

if (input_flag>1.5)

    string1=' tdoa_update_m1: after ';

end

figure, 
imagesc(abs(Z_ol));
title([string1, 'Z_ol']);

        figure,
        subplot(2,1,1)
        imagesc(abs(W_fom(:,:,1)));
        title([string1, 'W_fom']);
                        
        subplot(2,1,2)
        imagesc(abs(W_fom(:,:,2)));
        

    
end
