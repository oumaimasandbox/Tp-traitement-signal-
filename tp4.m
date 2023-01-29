
%% 

clear all
close all
clc
Te = 5*1e-4;
f1 = 500;
f2 = 400;
f3 = 50;
t = 0:Te:5-Te;
fe = 1/Te;
N = length(t);

fshift = (-N/2:N/2-1)*(fe/N);
f = (0:N-1)*(fe/N);
%notre signal
x = sin(2*pi*f1*t)+sin(2*f2*pi*t)+sin(2*pi*f3*t);
%on applique la transforme de fourrier
y = fft(x);


subplot(2,1,1)
plot(t,x)
xlim([0,0.5])

subplot(2,1,2)
plot(fshift, fftshift(abs(y)));

%% 

k = 1;

w=2*pi*f;
%on choisit 3 frequences de coupures
wc = 50;
wc1 = 500;
wc2 = 1000;
%la transmitance complexe 
h = (k*1j*((w)/wc))./(1+1j*((w)/wc));
h1 = (k*1j*((w)/wc1))./(1+1j*((w)/wc1));
h2 = (k*1j*((w)/wc2))./(1+1j*((w)/wc2));
%diagramme de bode en fct du gain
G = 20*log(abs(h));
G1 = 20*log(abs(h1));
G2 = 20*log(abs(h2));
%diagramme de bode en fct de la phase 
P = angle(h);
P1 = angle(h1);
P2 = angle(h2);

% subplot(3,1,1)
% semilogx(abs(h))
% % plot(abs(h))
% legend("Module de h(t)")
% 
% subplot(3,1,2)
% semilogx(f,G,f,G1,f,G2);
% title("Diagramme de Bode")
% xlabel("rad/s")
% ylabel("decibel")
% legend("G : wc=50","G1 : wc=500","G2 : wc=1000")
% 
% subplot(3,1,3)
% semilogx(f,P,f,P1,f,P2)
% legend("P","P1","P2")

%% 



Filtre1=h.*y;
Filtre2=h1.*y;
Filtre3=h2.*y;
%Retour au domaine temporel
x1=ifft(Filtre1,"symmetric");
x2=ifft(Filtre2,"symmetric");
x3=ifft(Filtre3,"symmetric");





%representation des signaux

subplot(4,1,1);
plot(t,x);
legend("signal initial")
xlim([0,0.25])
subplot(4,1,2);
plot(t,Filtre1);
legend("bruit wc=50")
xlim([0,0.25])
subplot(4,1,3);
plot(t,Filtre2);
legend("bruit wc=500")
xlim([0,0.25])
subplot(4,1,4);
plot(t,Filtre3);

legend("bruit wc=1000")
xlim([0,0.25])

%% 

% % debruitage
[data,fs]=audioread("test.wav");
Fe=11025;
Te=1/Fe;
N=length(data);
% t = 0:Te:(N-1)*Te;
plot(t,data)
grid on
% % 

%% 
[data,fs]=audioread("test.wav");
Fe=11025;
Te=1/Fe;
N=length(data);
t = 0:Te:5-Te;
f = (0:N-1)*(fs/N);
k = 1;
fc = 4700;%on a pas choisi 5000 parceque l'attenuation est de 0.7
%la transmitance complexe 
h =k./(1+1j*(f/fc).^2000;
%créer un filtre symétrique qui est nécessaire pour filtrer le signal 
% de manière symétrique et éviter des distorsions de phase dans le signal filtré.
h_filter = [h(1:floor(N/2)),flip(h(1:floor(N/2)))];

semilogx(f(1:floor(N/2)),abs( h(1:floor(N/2))),'linewidth',1.5)

y_filtr = y_trans(1:end-1).*h_filter;
sig_filtred= ifft(y_filtr,"symmetric");

plot(fshift(1:end-1),fftshift(abs(fft(sig_filtred))))
%plot(t,music-y_filtr )
% sound(sig_filtred,fs)

