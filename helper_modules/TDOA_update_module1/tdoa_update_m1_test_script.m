clear all
close all
clc

F_static=10; N_static=6; M_static=2; O_static=5; L_static=3; K_static=4;

Xtilde_fnm=rand(F_static, N_static, M_static)+1i*rand(F_static, N_static, M_static);
Z_ol=rand(O_static, L_static);
Y_lk=rand(L_static, K_static);
T_fk=rand(F_static, K_static);
V_nk=rand(N_static, K_static);
W_fom=rand(F_static, O_static, M_static);
expj_Phi_W_fom=rand(F_static, O_static, M_static)+1i*rand(F_static, O_static, M_static);
expj_Phi_S_fkn=rand(F_static, K_static, N_static)+1i*rand(F_static, K_static, N_static);

[Z_ol_out, W_fom_out] = tdoa_update_m1(Xtilde_fnm, Z_ol, Y_lk, T_fk, V_nk, W_fom, expj_Phi_W_fom, expj_Phi_S_fkn);