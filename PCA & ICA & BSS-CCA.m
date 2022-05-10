clc;
clear all;
close all;
%addpath E:\uni\arshad\term 2\EEG\HW_comp\2\Ex4;
load('Ex4.mat')

%%
for i=1:32
    figure(1)
    subplot(8,4,i)
    plot(X_org(i,:))
end

ss = sum(sum(X_org.^2));
[r,c] = size(X_org);
snr = [-5 -10 -15 -20];

n1 = sum(sum(X_noise_1.^2));
n2 = sum(sum(X_noise_2.^2));
n3 = sum(sum(X_noise_3.^2));
n4 = sum(sum(X_noise_4.^2));
n5 = sum(sum(X_noise_5.^2));

An1 = ((ss/n1)*(10.^(-snr./10))).^0.5;
An2 = ((ss/n2)*(10.^(-snr./10))).^0.5;
An3 = ((ss/n3)*(10.^(-snr./10))).^0.5;
An4 = ((ss/n4)*(10.^(-snr./10))).^0.5;
An5 = ((ss/n5)*(10.^(-snr./10))).^0.5;

%% SNR -5dB
X = X_org + An1(1)*X_noise_1 + An2(1)*X_noise_2 + An3(1)*X_noise_3 + An4(1)*X_noise_4 + An5(1)*X_noise_5;

% PCA
cov_X = cov(transpose(X));
[U,L]=eig(cov_X);
D = (L^-0.5)*transpose(U);
S_pca_5dB = D*X;
for i=1:16
    figure(2)
    subplot(8,2,i)
    plot(S_pca_5dB(i,:))
    figure(3)
    subplot(8,2,i)
    plot(S_pca_5dB(i+16,:))
end
A = inv(D);
SelSources = [1 10 11 12 13 14 18 20 28];
X_den_pca = A(:,SelSources) * S_pca_5dB(SelSources,:);
% for i=1:32
%     figure(4)
%     subplot(8,4,i)
%     plot(X_den_pca(i,:))
% end
%channel 13
figure(5)
subplot(3,1,1)
plot(X_den_pca(13,:))
title('X-denoise-PCA channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(6)
subplot(3,1,1)
plot(X_den_pca(24,:))
title('X-denoise-PCA channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_pca_5dB = (sum(sum((X_den_pca-X_org).^2)) / sum(sum(X_org.^2)))^0.5;
%%
% BSS-CCA
[S_CCA_5dB,A_x_CCA_5dB,autocor_CCA_5dB]= EMGsorsep(X);
for i=1:16
    figure(7)
    subplot(8,2,i)
    plot(S_CCA_5dB(i,:))
    figure(8)
    subplot(8,2,i)
    plot(S_CCA_5dB(i+16,:))
end
SelSources = [1 2 13 17 19 21 23 25 26 27 28 30 32];
X_den_CCA = A_x_CCA_5dB(:,SelSources) * S_CCA_5dB(SelSources,:);

% for i=1:32
%     figure(9)
%     subplot(8,4,i)
%     plot(X_den_CCA(i,:))
% end

%channel 13
figure(10)
subplot(3,1,1)
plot(X_den_CCA(13,:))
title('X-denoise-CCA channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(11)
subplot(3,1,1)
plot(X_den_CCA(24,:))
title('X-denoise-CCA channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

X_den_CCA(:,10240) = mean(X_den_CCA,2);
RRMSE_CCA_5dB = (sum(sum((X_den_CCA-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%%
%ICA-COM2R
[F,W,K]=COM2R(X,32);
S_COM2R_5dB = W*X;
for i=1:16
    figure(12)
    subplot(8,2,i)
    plot(S_COM2R_5dB(i,:))
    figure(13)
    subplot(8,2,i)
    plot(S_COM2R_5dB(i+16,:))
end
A = inv(W);
SelSources = [1 3 4 5 7 9 10 11 17 19 21 28 30];
X_den_COM2R = A(:,SelSources) * S_COM2R_5dB(SelSources,:);
% for i=1:32
%     figure(14)
%     subplot(8,4,i)
%     plot(X_den_COM2R(i,:))
% end
%channel 13
figure(15)
subplot(3,1,1)
plot(X_den_COM2R(13,:))
title('X-denoise-COM2R channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(16)
subplot(3,1,1)
plot(X_den_COM2R(24,:))
title('X-denoise-COM2R channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_COM2R_5dB = (sum(sum((X_den_COM2R-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%%
%ICAsobi
[H,S_sobi_5dB,D]=sobi(X,32,32);
for i=1:16
    figure(17)
    subplot(8,2,i)
    plot(S_sobi_5dB(i,:))
    figure(18)
    subplot(8,2,i)
    plot(S_sobi_5dB(i+16,:))
end

SelSources = [2 4 6 17 23 32];
X_den_sobi = H(:,SelSources) * S_sobi_5dB(SelSources,:);

% for i=1:32
%     figure(19)
%     subplot(8,4,i)
%     plot(X_den_CCA(i,:))
% end

%channel 13
figure(20)
subplot(3,1,1)
plot(X_den_sobi(13,:))
title('X-denoise-sobi channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(21)
subplot(3,1,1)
plot(X_den_sobi(24,:))
title('X-denoise-sobi channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_sobi_5dB = (sum(sum((X_den_sobi-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%% SNR -20dB
X = X_org + An1(4)*X_noise_1 + An2(4)*X_noise_2 + An3(4)*X_noise_3 + An4(4)*X_noise_4 + An5(4)*X_noise_5;

% PCA
cov_X = cov(transpose(X));
[U,L]=eig(cov_X);
D = (L^-0.5)*transpose(U);
S_pca_20dB = D*X;
for i=1:16
    figure(22)
    subplot(8,2,i)
    plot(S_pca_20dB(i,:))
    figure(23)
    subplot(8,2,i)
    plot(S_pca_20dB(i+16,:))
end
A = inv(D);
SelSources = [1 5 10 11 27 29];
X_den_pca = A(:,SelSources) * S_pca_20dB(SelSources ,:);
% for i=1:32
%     figure(24)
%     subplot(8,4,i)
%     plot(X_den_pca(i,:))
% end
%channel 13
figure(25)
subplot(3,1,1)
plot(X_den_pca(13,:))
title('X-denoise-PCA channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')
%channel 24
figure(26)
subplot(3,1,1)
plot(X_den_pca(24,:))
title('X-denoise-PCA channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_pca_20dB = (sum(sum((X_den_pca-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%%
% BSS-CCA
[S_CCA_20dB,A_x_CCA_20dB,autocor_CCA_5dB]= EMGsorsep(X);
for i=1:16
    figure(27)
    subplot(8,2,i)
    plot(S_CCA_20dB(i,:))
    figure(28)
    subplot(8,2,i)
    plot(S_CCA_20dB(i+16,:))
end
SelSources = [1 3 12 19 20 23 25 27 29 30 31];
X_den_CCA = A_x_CCA_20dB(:,SelSources) * S_CCA_20dB(SelSources,:);

% for i=1:32
%     figure(29)
%     subplot(8,4,i)
%     plot(X_den_CCA(i,:))
% end

%channel 13
figure(30)
subplot(3,1,1)
plot(X_den_CCA(13,:))
title('X-denoise-CCA channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(31)
subplot(3,1,1)
plot(X_den_CCA(24,:))
title('X-denoise-CCA channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

X_den_CCA(:,10240) = mean(X_den_CCA,2);
RRMSE_CCA_20dB = (sum(sum((X_den_CCA-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%%
%ICA-COM2R
[F,W,K]=COM2R(X,32);
S_COM2R_20dB = W*X;
for i=1:16
    figure(32)
    subplot(8,2,i)
    plot(S_COM2R_20dB(i,:))
    figure(33)
    subplot(8,2,i)
    plot(S_COM2R_20dB(i+16,:))
end
A = inv(W);
SelSources = [2 4 5 7 8 10 17 20 21 24 31 32];
X_den_COM2R = A(:,SelSources) * S_COM2R_20dB(SelSources ,:);
% for i=1:32
%     figure(34)
%     subplot(8,4,i)
%     plot(X_den_COM2R(i,:))
% end
%channel 13
figure(35)
subplot(3,1,1)
plot(X_den_COM2R(13,:))
title('X-denoise-COM2R channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-rg channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(36)
subplot(3,1,1)
plot(X_den_COM2R(24,:))
title('X-denoise-COM2R channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_COM2R_20dB = (sum(sum((X_den_COM2R-X_org).^2)) / sum(sum(X_org.^2)))^0.5;

%%
%ICAsobi
[H,S_sobi_20dB,D]=sobi(X,32,32);
for i=1:16
    figure(37)
    subplot(8,2,i)
    plot(S_sobi_20dB(i,:))
    figure(38)
    subplot(8,2,i)
    plot(S_sobi_20dB(i+16,:))
end

SelSources = [6 17 18 19 32];
X_den_sobi = H(:,SelSources) * S_sobi_20dB(SelSources,:);

% for i=1:32
%     figure(39)
%     subplot(8,4,i)
%     plot(X_den_CCA(i,:))
% end

%channel 13
figure(40)
subplot(3,1,1)
plot(X_den_sobi(13,:))
title('X-denoise-sobi channel 13')
subplot(3,1,2)
plot(X_org(13,:))
title('X-org channel 13')
subplot(3,1,3)
plot(X(13,:))
title('X-noisy channel 13')

%channel 24
figure(41)
subplot(3,1,1)
plot(X_den_sobi(24,:))
title('X-denoise-sobi channel 24')
subplot(3,1,2)
plot(X_org(24,:))
title('X-org channel 24')
subplot(3,1,3)
plot(X(24,:))
title('X-noisy channel 24')

RRMSE_sobi_20dB = (sum(sum((X_den_sobi-X_org).^2)) / sum(sum(X_org.^2)))^0.5;
