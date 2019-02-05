function [Xhat_fnm, W_fom, T_fk, V_nk, Phi_S_fnk, S_fnl, X_R, Freqxr, Timexr, X_L, Freqxl, Timexl]=nguyen_2015_matlab()

clear all
close all
clc

diary('nguyen_2015_matlab_diary.txt');

load('input_LR_instruments_half_plane.mat');
%load('input_LR_instruments_no_reverb.mat');
%load('/home/tung/Documents/School/thesis2/mycode/sawada_may2013/code/input_LR.mat');

%load('/mnt/asus_share/School/thesis/mycode/matlab/sawada_may2013/code/input_LR.mat');
%load('/home/tung/Documents/School/thesis/mycode/matlab/sawada_may2013/code/input_LR.mat');

%input_left:101599x1 double
%input_right:101599x1 double

M=2; %channel #.

%L=9;
% L=30;
% K=60;
L=24;
K=24;

A=36;       %azimuth
E=8;        %elevation
% O=A*(E-1)+2;      %total # look directions: 254;
%O=110;
O=55;

%%STFTs----------------------------------------------------------------------------------

curr_dir=pwd; %save the current directory as a string

%cd /mnt/asus_share/School/thesis/mycode/matlab/sawada_may2013/code/mex_eu_nmf/istft
cd /home/tung/Documents/School/thesis/mycode/matlab/sawada_may2013/code/mex_eu_nmf/istft

%pwd
wlen=2*512;
% wlen=2*1024;
h=wlen/4;
nfft=wlen;
Fs=16000;


[X_R, Freqxr, Timexr]=spectrogram(input_right,hamming(1024),256,512,16000);
display(length(Freqxr));
[X_R, Freqxr, Timexr]=stft(input_right, wlen, h, nfft, Fs);

[X_L, Freqxl, Timexl]=spectrogram(input_left,hamming(1024),256,512,16000);
display(length(Timexr));
[X_L, Freqxl, Timexl]=stft(input_left, wlen, h, nfft, Fs);

%spectrograms of Y1, Y2, Y3
[Y1_stft_FxN, Freq_Y1, Time_Y1]=stft(Y1, wlen, h, nfft, Fs);

%figure, 
%imagesc(abs(Y1_stft_FxN));

%figure,
%imagesc(angle(Y1_stft_FxN));

[Y2_stft_FxN, Freq_Y2, Time_Y2]=stft(Y2, wlen, h, nfft, Fs);

%figure, 
%imagesc(abs(Y2_stft_FxN));

%figure,
%imagesc(angle(Y2_stft_FxN));

[Y3_stft_FxN, Freq_Y3, Time_Y3]=stft(Y3, wlen, h, nfft, Fs);

%figure, 
%imagesc(abs(Y3_stft_FxN));

%figure,
%imagesc(angle(Y3_stft_FxN));

cd(curr_dir) %cd back to the original directory

%end
%STFTs-----------------------------------------------------------------------------------

%%Randomize Matrices & initialize dimensions

F=length(Freqxr);
display(F);

N=length(Timexr);
display(N);

Xtilde_fnm1=zeros(F,N,M);

Xtilde_fnm1(:,:,1)=X_R;
Xtilde_fnm1(:,:,2)=X_L;

%Preprocess by taking the amplitude square rooted of Xtilde_fnm1:
for f=1:F
    
   for n=1:N
    
       for m=1:M
       
          Xtilde_fnm1(f,n,m)=sqrt(abs(Xtilde_fnm1(f,n,m)))*exp(1i*angle(Xtilde_fnm1(f,n,m))); 
           
       end
   end
    
end

%H_flm=rand(F,L,M);
% T_fk=rand(F,K); %+ones(F,K);
T_fk=sqrt(1/K)*ones(F,K)+(1/64)*randn(F,K);

%Scale T_fk by the l1-norm for k=1:K

%populate V_nk
% V_nk=rand(N,K); %+ones(N,K);
V_nk=sqrt(1/K)*ones(N,K)+(1/64)*randn(N,K);

%Phi_H_flm=2*pi*(rand(F,L,M)-0.5);

%Phi_S_fnk=2*pi*(rand(F,N,K)-0.5);

Phi_S_nkf=2*pi*(rand(N,K,F)-0.5);
Phi_S_fkn=zeros(F,K,N);

for n=1:N
   
    Phi_S_fkn(:,:,n)=transpose(squeeze(Phi_S_nkf(n,:,:)));
    
end

expj_S_nkf=zeros(N,K,F);
expj_S_fkn=zeros(F,K,N);


expj_S_nkf=exp(1i*Phi_S_nkf);
expj_S_fkn=exp(1i*Phi_S_fkn);
%BUGFIX: populating these tensors using this intendedly overloaded code
%doesn't seem to work well for large-sized tensors. populate it with a
%triple for loop instead.

% for f=1:F
%    
%     for k=1:K
%        
%         for n=1:N
%         
%             expj_S_nkf(n,k,f)=exp(1i*Phi_S_nkf(n,k,f));
%             expj_S_fkn(f,k,n)=exp(1i*Phi_S_fkn(f,k,n));
%             
%         end
%         
%     end
%     
% end

%save('Phi_S_fnkl.mat', 'Phi_S_fnkl');

%S_fnl=zeros(F,N,L);

% W_fom=rand(F,O,M);
W_fom=sqrt(1/O)*ones(F,O,M)+(1/64)*randn(F,O,M);

%Normalize each Z_ol column to 1/O to start. 
Z_ol=sqrt(1/O)*ones(O, L)+(1/64)*randn(O, L); 
% Z_ol=rand(O,L);

%Compute its new mean..

% Z_ol_new_mean=(sum(sum(Z_ol))/numel(Z_ol));

%compute its current mean
% W_fom_current_mean=sum(sum(sum(W_fom)))/numel(W_fom);
% 
% W_fom=W_fom-W_fom_current_mean*ones(F,O,M);    %shift the mean to become a zero mean tensor
% 
% W_fom=W_fom+(1/(Z_ol_new_mean*O))*ones(F,O,M);

%Normalize each Y_lk column to 1/L to start. 
Y_lk=(1/L)*ones(L,K)+(1/64)*randn(L,K);

v_sound=340.3; 
mic_coordinates=zeros(3, M); 
k0_dummy=zeros(3,1);

%Set the mic coordinate vectors
mic_coordinates(:, 1)=[ 0; 0.5; 0];
mic_coordinates(:, 2)=[ 0; -0.5; 0];

Phi_W_fom=zeros(F,O,M);

N_stft=nfft; 

% Uncomment if you want a new set of random variables. For example, if you
% changed the length of O. 
u_rand=rand(O,1);
v_rand=rand(O,1);

save('u_v_rand.mat', 'u_rand', 'v_rand');

load('u_v_rand.mat');

for o=1:O

    %elevation of matlab and this tutorial are defined differently!: http://www.mathworks.com/help/phased/ug/spherical-coordinates.html#bs98aup-1
    %http://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
    %will switch the tutorial's way of populating the cartesian coordinates
    phi=acos(2*v_rand(o)-1);
    %theta=2*pi*u_rand(o);
    theta=(u_rand(o)-0.5)*(pi); %a uniform distribution on the interval -0.5 to 0.5 then multiplied by pi should have max value pi*0.5, ie: pi/2
    
    k0_dummy(1)=sin(phi)*cos(theta);    %expect that the output x coordinate be strictly positive. This involves I think restricting the azimuth (theta) to be + or - 90deg (pi/2) instead of + or - 180deg (pi)
    k0_dummy(2)=sin(phi)*sin(theta);
    k0_dummy(3)=cos(phi);    
    
    %matlab way: wrong, assuming random initalizations of phi and theta as
    %a function of u and v....
%     k0_dummy(1)=cos(theta)*cos(phi);
%     k0_dummy(2)=cos(theta)*sin(phi);
%     k0_dummy(3)=sin(theta);        
    
    for m=1:M
    
        Tau_k0=-sum(k0_dummy.*squeeze(mic_coordinates(:, m)));        
        
        for f=1:F
            
            Phi_W_fom(f, o, m)=(2*pi*(f-1)*Fs*Tau_k0)/(N_stft*F);            
            
        end
        
    end
    
end

expj_W_fom=zeros(F,O,M);        
expj_W_fom=exp(1i*Phi_W_fom);        
%BUGFIX::: ^exp(.) overloading notation doesn't seem to work properly for large sized
%matrices.

%Therefore, populate expj_W_fom with a for loop
% for m=1:M
%    
%     for o=1:O
%        
%         for f=1:F
%         
%             expj_W_fom(f,o,m)=exp(1i*Phi_W_fom(f,o,m));
%             
%         end
%     end
%     
% end
        
%decrease O momentarily for debugging

% O2=10;
% Phi_W_fom=rand(F,O2,M);        
% W_fom=rand(F,O2,M);  
% Z_ol=(1/O2)*ones(O2, L); 

%Generate the look direction magnitudes and phases for o=1:O, stringing the matrix (azimuth, elevation indices) out as
%a vector/sequence. 

% figure, plot_T_fk(T_fk)
%figure, hist(T_fk)

% figure, plot_V_nk(V_nk)
%figure, hist(V_nk)

% figure, plot_W_fom(W_fom)
%figure, hist(W_fom)

% figure, plot_Z_ol(Z_ol)
%figure, hist(Z_ol)

% figure, plot_Y_lk(Y_lk)
%figure, hist(Y_lk)

display(M)

% figure, imagesc(squeeze(Phi_W_fom(:,:,1)));
% 
% figure, imagesc(squeeze(Phi_W_fom(:,:,2)));

%load_and_plot_converged_stuff();

[Xhat_fnm, W_fom, Z_ol, Y_lk, T_fk, V_nk, expj_S_nkf, expj_S_fkn]=nguyen_2015_mex(Xtilde_fnm1, 0, W_fom, Z_ol, Y_lk, T_fk, V_nk, Phi_W_fom, Phi_S_nkf, M, F, N, K, L, O, Phi_S_fkn, expj_S_nkf, expj_S_fkn, expj_W_fom);
%[Xhat_fnm, W_fom, Z_ol, Y_lk, T_fk, V_nk, Phi_S_fnk]=nguyen_2015_mex(Xtilde_fnm1, 0, W_fom, Z_ol, Y_lk, T_fk, V_nk, Phi_W_fom, Phi_S_fnk, M, F, N, K, L, 5);

% close all
% clc 
% 
% addpath('/home/tung/Documents/School/thesis2/mycode/matlab/nguyen_2015_2/code_vectorised_16/bottom_up/mex')
% 
% Phi_S_nkf=angle(expj_S_nkf);
% 
% [Xhat_fnm, W_fom, Z_ol, Y_lk, T_fk, V_nk, expj_S_nkf, expj_S_fkn]=nguyen_2015_bottom_up_mex(Xtilde_fnm1, 0, W_fom, Z_ol, Y_lk, T_fk, V_nk, Phi_W_fom, Phi_S_nkf, M, F, N, K, L, O, Phi_S_fkn, expj_S_nkf, expj_S_fkn, expj_W_fom);

save('nguyen_2015_output_file_instruments_half_plane.mat', 'Xhat_fnm', 'W_fom', 'Z_ol', 'Y_lk', 'T_fk', 'V_nk', 'expj_S_nkf', 'expj_S_fkn', 'M', 'F', 'N', 'K', 'L', 'O', 'expj_W_fom', 'Xtilde_fnm1');

end

function load_and_plot_converged_stuff()

load('nguyen_2015_output_file_instruments_half_plane');

plot_T_fk(T_fk, 0);

plot_V_nk(V_nk);

plot_W_fom(W_fom);

plot_Z_ol(Z_ol);

plot_Y_lk(Y_lk);

plot_Phi_S_fnk(expj_S_fkn, 24, 1);

end