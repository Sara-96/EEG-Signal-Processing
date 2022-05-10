clc;
clear all;
close all;
load('Electrodes.mat')
load('NewData1.mat')
EEG_Sig1 = EEG_Sig;
load('NewData2.mat')
EEG_Sig2 = EEG_Sig;
load('NewData3.mat')
EEG_Sig3 = EEG_Sig;
load('NewData4.mat')
EEG_Sig4 = EEG_Sig;
%%
X = EEG_Sig1;
plotEEG
title('Sig1')

X = EEG_Sig2;
plotEEG
title('Sig2')

X = EEG_Sig3;
plotEEG
title('Sig3')

X = EEG_Sig4;
plotEEG
title('Sig4')

load('Electrodes') ;
[elocsX,elocsY] = pol2cart(pi/180*[Electrodes.theta],[Electrodes.radius]);

%% Sig1
fs=250;
t=0:1/fs:1250/fs;
[F1,W1,K1]=COM2R(EEG_Sig1,21);
S1 = W1*EEG_Sig1;
for i=1:21
    figure(5)
    subplot(7,3,i)
    plot(t,S1(i,:))
    title(['time / channel ',num2str(i)])
    
    [pxx,f] = pwelch(S1(i,:),fs);
    figure(6)
    subplot(7,3,i)
    plot(f*fs/pi/2,10*log10(pxx))
    title(['frequency / channel ',num2str(i)])

    figure(7)
    subplot(7,3,i)
    plottopomap(elocsX, elocsY, Electrodes.labels, F1(:,i));
    title(['channel ',num2str(i)])
end

A = inv(W1);
SelSources = [10];
X_denoise_1 = A(:,SelSources)*S1(SelSources,:);
X = X_denoise_1;
plotEEG
title('Sig1-denoise')
%% Sig2
fs=250;
t=0:1/fs:2500/fs;
[F2,W2,K2]=COM2R(EEG_Sig2,21);
S2 = W2*EEG_Sig2;
for i=1:21
    figure(9)
    subplot(7,3,i)
    plot(t,S2(i,:))
    title(['time / channel ',num2str(i)])
    
    [pxx,f] = pwelch(S2(i,:),fs);
    figure(10)
    subplot(7,3,i)
    plot(f*fs/pi/2,10*log10(pxx))
    title(['frequency / channel ',num2str(i)])

    figure(11)
    subplot(7,3,i)
    plottopomap(elocsX, elocsY, Electrodes.labels, F2(:,i));
    title(['channel ',num2str(i)])
end

A = inv(W2);
SelSources = [3 17];
X_denoise_2 = A(:,SelSources)*S2(SelSources,:);
X = X_denoise_2;
plotEEG
title('Sig2-denoise')

%% Sig3
fs=250;
t=0:1/fs:1250/fs;
[F3,W3,K3]=COM2R(EEG_Sig3,21);
S3 = W3*EEG_Sig3;
for i=1:21
    figure(13)
    subplot(7,3,i)
    plot(t,S3(i,:))
    title(['time / channel ',num2str(i)])
    
    [pxx,f] = pwelch(S3(i,:),fs);
    figure(14)
    subplot(7,3,i)
    plot(f*fs/pi/2,10*log10(pxx))
    title(['frequency / channel ',num2str(i)])

    figure(15)
    subplot(7,3,i)
    plottopomap(elocsX, elocsY, Electrodes.labels, F3(:,i));
    title(['channel ',num2str(i)])
end
A = inv(W3);
SelSources = [6 9 13 14 15];
X_denoise_3 = A(:,SelSources)*S3(SelSources,:);
X = X_denoise_3;
plotEEG
title('Sig3-denoise')
%% Sig4
fs=250;
t=0:1/fs:1250/fs;
[F4,W4,K4]=COM2R(EEG_Sig4,21);
S4 = W4*EEG_Sig4;
for i=1:21
    figure(17)
    subplot(7,3,i)
    plot(t,S4(i,:))
    title(['time / channel ',num2str(i)])
    [pxx,f] = pwelch(S4(i,:),fs);
    figure(18)
    subplot(7,3,i)
    plot(f*fs/pi/2,10*log10(pxx))
    title(['frequency / channel ',num2str(i)])

    figure(19)
    subplot(7,3,i)
    plottopomap(elocsX, elocsY, Electrodes.labels, F4(:,i));
    title(['channel ',num2str(i)])
end
A = inv(W4);
SelSources = [13 14 18 19 20 21];
X_denoise_4 = A(:,SelSources)*S4(SelSources,:);
X = X_denoise_4;
plotEEG
title('Sig4-denoise')