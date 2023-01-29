clear all 
close all 
clc


load('ecg.mat')
Fe=500;
Te=1/Fe;
N=length(ecg);
t = 0:Te:(N-1)*Te;
subplot(3,2,1)
plot(t,ecg)
grid on
title(" représentation graphique de l’activation électrique du cœur")
xlabel("t")
ylabel("ECG")
xlim([0.5 1.5]);

%%
f=(0:N-1)*(Fe/N);
fshift=(-N/2:N/2-1)*(Fe/N);
y = fft(ecg);
% subplot(3,2,2)
plot(f,fftshift(abs(y)));
grid on
title(" représentation graphique de la transformée de fourier du signal ECG")
xlabel("f")
ylabel("Tfd")
%%
% conception du filtre
pass_haut=ones(size(ecg));
fc=0.5;
index_fc= ceil((fc*N)/Fe);
pass_haut(1:index_fc)= 0;
pass_haut(N-index_fc+1:N) = 0;
 
plot(f,pass_haut,"linewidth",1.5)
%% Filtrage

ecg1_freq = pass_haut.*y; 
ecg1 = ifft(ecg1_freq,"symmetric");

%% plot filterd signal 
subplot(3,2,2)
plot(t,ecg1)
hold on
plot(t,ecg+3)
hold on 
plot(t,ecg-ecg1+1.5)



%% Filtrage fréquence 50%

pass_notch=ones(size(ecg));
%crée un vecteur  qui a la même taille que le signal ECG et est rempli de 1.
fc2=50;
index_fc2= ceil((fc2*N)/Fe)+1;
%calcule l'index  de fréquence 50 Hz dans
% la représentation de domaine fréquentiel du signal ECG en utilisant la fréquence d'échantillonnage
pass_notch(index_fc2)= 0;
pass_notch(N-index_fc2+1) = 0;

%% 


ecg2_freq = pass_notch.*fft(ecg1); 
%on recupere le signal dand le domaine temporel
ecg2 = ifft(ecg2_freq,"symmetric");


subplot(3,1,1)

plot(t,ecg1,'r')
 grid on
subplot(3,1,2)
plot(t,ecg1-ecg2)
 grid on 
subplot(3,1,3)

plot(t,ecg2,'b')
grid on


%% Filtrage pass_bas 

pass_bas=zeros(size(ecg));
fc3 = 20;
index_fc3 = ceil((fc3*N)/Fe);
pass_bas(1:index_fc3)= 1;
pass_bas(N-index_fc3+1:N) = 1;


%% Filtrage 2 

ecg3_freq = pass_bas.*fft(ecg2); 
ecg3 = ifft(ecg3_freq,"symmetric");

subplot(313)
plot(t,ecg3)
grid on
subplot(311)
plot(t,ecg2)
grid on
subplot(312)
plot(t,ecg2-ecg3)
grid on
 
% plot(t,ecg-ecg3)
% 
% 
% [c,lags] = xcorr(ecg3,ecg3)

%% 
% Charger le signal ECG traité 
ecg_signal = load('ecg.mat');

% Définir l'intervalle de recherche pour la fréquence cardiaque
min_bc = 40; 
max_bc = 220; 

% Calculer l'autocorrélation du signal ECG
[acf,lags] = xcorr(ecg_signal,ecg_signal);

% Trouver la fréquence cardiaque en se basant sur l'autocorrélation
[max_corr, max_index] = max(acf);
heart_rate = 60/(lags(max_index));

% Vérifier si la fréquence cardiaque est dans l'intervalle de recherche
if heart_rate > min_bc && heart_rate < max_bc
    disp(['Fréquence cardiaque : ', num2str(heart_rate), ' battements par minute']);
else
    disp('Fréquence cardiaque non détectée');
end




