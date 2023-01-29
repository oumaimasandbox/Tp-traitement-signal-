clear all
close all
clc

fe = 1e4;
Te = 1/fe;
N = 5000; 
% generer un vecteur t qui contient des valeurs discretes et on applique le
% cosinus sur le vecteur t
  t = 0: Te: (N -1)*Te;
 %un signal enchantillone avec une frequence de 10^4
x = 1.2*cos(2*pi*440*t + 1.2) + 3*cos(2*pi*550*t) + 0.6*cos(2*pi*2500*t);
subplot(3, 2, 1)
plot(t,x,".")
% % Application de la tfd
%on configure l'axe frequentiel avec un pas de descretisation de fe/N
fshift = 0: fe/N : (N -1)*(fe/N);
x = 1.2*cos(2*pi*440*t + 1.2) + 3*cos(2*pi*550*t) + 0.6*cos(2*pi*2500*t);
subplot(3, 2, 2)
%on represente le module du spectre d'amplitude
y= abs(fft(x));
plot(fshift,y)

%on va injecter le bruit dans le signal
subplot(3, 2, 3)
xbruit = x + 2*randn(size(x));
plot(t,2*randn(size(x)))

subplot(3, 2, 4)
ybruit = abs(fft(xbruit));
plot(fshift,fftshift(ybruit))
%on va intensifier le bruit dans le signal en augmentant le coefficients

subplot(3, 2, 5)
xbruit2 = x + 50*randn(size(x));
plot(t,50*randn(size(x)))

subplot(3, 2, 6)
ybruit3 = abs(fft(xbruit2));
plot(fshift,fftshift(ybruit3))

%% clear all
close all
clc 

fe=10000;
Te=1/fe;
N=5000;

t = 0:Te:(N-1)*Te;
x = 1.2*cos(2*pi*440*t+1.2)+3*cos(2*pi*550*t)+0.6*cos(2*pi*2500*t);
x2 = 1.2*cos(2*pi*440*t+1.2)+3*cos(2*pi*550*t);

f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*(fe/N);
y = fft(x);

% Conception du filtre 
  %initialisation du filtre
pass_bas = zeros(size(x));
 %definir la frequence de coupure
fc = 2000;
%l'index de la frequence de coupure

index_fc = ceil((fc*N)/fe);
%ceil Arrondit chaque élément de X à l'entier le plus proche, supérieur ou égal à cet élément.
pass_bas(1:index_fc) = 1;
pass_bas(N-index_fc+1:N) = 1;
% plot(f,pass_bas,"linewidth",1.5)
legend("filtre")
% Filtrage

x_filter_freq = pass_bas.*y;
x_filter_temp = ifft(x_filter_freq,"symmetric");

%  plot(fshift, fftshift(abs(fft(x_filter_temp))));

 plot(t,x2-x_filter_temp)








%% Analyse fréquentielle du chant du rorqual bleu
%on enregistre l'audiodans une variable nommée "data" et la fréquence d'échantillonnage dans une variable nommée "fs"
[data,fs]=audioread("bluewhale (1).wav");
sound(data,fs);

son1 = data(2.45e4: 3.10e4);
%sound(data,fs);

plot(son1);
xlabel('t');
ylabel('son');
title('Représentation temporelle du son roqual bleu')
grid on

%% 
    taille=size(data);
  %On applique la transformation de fourier rapide sur le chant 
     Schant = fft(data);
    %Densité spectrale du Chant
    Densite_spectrale_chant = abs(Schant).^2/taille;
    f = (0:floor(taille/2))*(fs/taille)/10;
    plot(f,Densite_spectrale_chant(1:floor(taille/2)+1));
   
    %% 

%question 3:
N=length(son);
x = abs(fft(son)).^2/N; 
f = (0:floor(N/2))*(Fe/N)/10;
plot(f,x(1:floor(N/2)+1));
title('Le signal densité spectrale de puissance du signal')
grid on
