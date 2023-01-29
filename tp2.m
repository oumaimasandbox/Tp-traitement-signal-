clear all
close all
clc
%on enregistre l'audiodans une variable nommée "data" et la fréquence d'échantillonnage dans une variable nommée "fs"
[data,fs]=audioread("Centre.mp3");
%On trace le signal
Te=1/fs;
N=length(data);
t = 0:Te:(N-1)*Te;
plot(t,data)
grid on
%pour ecouter l'audio 
%   sound(data,fs)
%On modifie la frequence d'echantillonage pour acceler/ ou ralentir l'audio

audio_ralentie = resample(data,fs*4,fs);
% sound(y2,fs);
audio_accelere= resample(data,fs/2,fs);
% sound(y2,fs/2);
%% diviser l'enregistrement puis le restiteuer

 
stem(data)
%pour obsever du signal enregistré x qui correspondent à chaque morceau. 


riennesertde=data(40992:98000)
%on segmente le signal 
% stem(riennesertde)
% sound(riennesertde,fs)
 
 courir=data(98000:140000);
%  sound(courir,fs)

 ilfaut=data(140000:160999);
%   sound(ilfaut,fs)
 
 partirapoint=data(160999:233384);
%  sound(partirapoint,fs)
 
parole=[riennesertde;courir;ilfaut;partirapoint];
sound(parole,fs)


%% Synthèse et analyse spectrale dune gamme de musique

 

 

fe=8192;
te=1/fe;
N=5000;
t=(0:N-1)*te;
do=10*cos(2*pi*262*t);
%sound(do,fe)
re=10*cos(2*pi*294*t);
%sound(re,fe)
mi=10*cos(2*pi*330*t);
%sound(mi,fe)
fa=10*cos(2*pi*349*t);
%sound(re,fa)
sol=10*cos(2*pi*392*t);
%sound(sol,fe)
la=10*cos(2*pi*440*t);
%sound(la,fe)
si=10*cos(2*pi*494*t);
%sound(si,fe)
do2=10*cos(2*pi*523*t);
%sound(do2,fe)
musique=[do,re,mi,fa,sol,la,si,do2];
 sound(musique,fe)
 
%% Spectre de la gamme de musique





f=(0:N-1)*(fs/N);
signalAnalyzer(musique);


spectre_musique=fft(musique);

%% Approximation du spectre d’un signal sinusoïdal à temps continu par FFT

 subplot(2,1,1)
    plot(fshift,fftshift(abs(y)));
    legend("Represenation du spectre d'une Octave");
    xlabel("f");
    ylabel("A");
    subplot(2,1,2)
    sig = 20*log(fftshift(abs(y)));

    plot(fshift,sig);
    legend("Represenation du spectre d'une Octave en dB");
    xlabel("f");
    ylabel("A");

