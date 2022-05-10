% HW3 - 97200216
clc; clear all; close all;
%%
load('HW3_Ex2.mat')
fs = 100;
N = length(X_org);
t = (0:(N-1))/fs;

%%
T_on = zeros(1,N);
for i=1:32
    a = find(abs(X_org(i,:))>3);
    T_on(a) = 1;
end

ss = sum(sum(X_org.^2));
snr = [-5 -10 -20];

n1 = sum(sum(X_noise_1.^2));
n2 = sum(sum(X_noise_2.^2));

An1 = ((ss/n1)*(10.^(-snr./10))).^0.5;
An2 = ((ss/n2)*(10.^(-snr./10))).^0.5;

%% SNR -5dB Noise1
X_noisy = X_org + An1(1)*X_noise_1;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_5dB_noise1 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_5dB_GEVD_noise1 = (sum(sum((X_org-X_den_5dB_noise1).^2))/sum(sum(X_org.^2)))^0.5;

figure(1)
subplot(3,1,1)
plot(t,X_den_5dB_noise1(13,:))
title('denoised channel 13 - 5dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 5dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 5dB')

figure(2)
subplot(3,1,1)
plot(t,X_den_5dB_noise1(24,:))
title('denoised channel 24 - 5dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 5dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 5dB')

%% SNR -10dB
X_noisy = X_org + An1(2)*X_noise_1;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_10dB_noise1 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_10dB_GEVD_noise1 = (sum(sum((X_org-X_den_10dB_noise1).^2))/sum(sum(X_org.^2)))^0.5;

figure(3)
subplot(3,1,1)
plot(t,X_den_10dB_noise1(13,:))
title('denoised channel 13 - 10dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 10dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 10dB')

figure(4)
subplot(3,1,1)
plot(t,X_den_10dB_noise1(24,:))
title('denoised channel 24 - 10dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 10dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 10dB')
%% SNR -20dB
X_noisy = X_org + An1(3)*X_noise_1;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_20dB_noise1 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_20dB_GEVD_noise1 = (sum(sum((X_org-X_den_20dB_noise1).^2))/sum(sum(X_org.^2)))^0.5;

figure(5)
subplot(3,1,1)
plot(t,X_den_20dB_noise1(13,:))
title('denoised channel 13 - 20dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 20dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 20dB')

figure(6)
subplot(3,1,1)
plot(t,X_den_20dB_noise1(24,:))
title('denoised channel 24 - 20dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 20dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 20dB')

%% SNR -5dB Noise2
X_noisy = X_org + An2(1)*X_noise_2;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_5dB_noise2 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_5dB_GEVD_noise2 = (sum(sum((X_org-X_den_5dB_noise2).^2))/sum(sum(X_org.^2)))^0.5;

figure(7)
subplot(3,1,1)
plot(t,X_den_5dB_noise2(13,:))
title('denoised channel 13 - 5dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 5dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 5dB')

figure(8)
subplot(3,1,1)
plot(t,X_den_5dB_noise2(24,:))
title('denoised channel 24 - 5dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 5dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 5dB')

%% SNR -10dB
X_noisy = X_org + An2(2)*X_noise_2;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_10dB_noise2 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_10dB_GEVD_noise2 = (sum(sum((X_org-X_den_10dB_noise2).^2))/sum(sum(X_org.^2)))^0.5;

figure(9)
subplot(3,1,1)
plot(t,X_den_10dB_noise2(13,:))
title('denoised channel 13 - 10dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 10dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 10dB')

figure(10)
subplot(3,1,1)
plot(t,X_den_10dB_noise2(24,:))
title('denoised channel 24 - 10dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 10dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 10dB')

% for i=1:16
%     figure(1)
%     subplot(8,2,i)
%     plot(t,X_org(i,:),t,X_den(i,:))
%     legend('org','den')
%     
%     figure(2)
%     subplot(8,2,i)
%     plot(t,X_org(i+16,:),t,X_den(i+16,:))
%     legend('org','den')
% end
%% SNR -20dB
X_noisy = X_org + An2(3)*X_noise_2;
Cx_org = cov(X_noisy');
a = zeros(32,32);
m = 0;
for i=1:N
    if(T_on(i)==1)
        a = a + (X_noisy(:,i)-mean(X_noisy,2))*transpose(X_noisy(:,i)-mean(X_noisy,2));
        m = m+1;
    end
end
Px_on = a ./ m;
[V_on,D_on] = eig(Px_on,Cx_org);
S_on = V_on' * X_noisy;
S_denoised_on = S_on;
S_denoised_on(1:31,:) = 0;
X_den_20dB_noise2 = inv(V_on') * S_denoised_on + repmat(mean(X_noisy,2),1,N);
error_20dB_GEVD_noise2 = (sum(sum((X_org-X_den_20dB_noise2).^2))/sum(sum(X_org.^2)))^0.5;

figure(11)
subplot(3,1,1)
plot(t,X_den_20dB_noise2(13,:))
title('denoised channel 13 - 20dB')
subplot(3,1,2)
plot(t,X_org(13,:))
title('original channel 13 - 20dB')
subplot(3,1,3)
plot(t,X_noisy(13,:))
title('noisy channel 13 - 20dB')

figure(12)
subplot(3,1,1)
plot(t,X_den_20dB_noise2(24,:))
title('denoised channel 24 - 20dB')
subplot(3,1,2)
plot(t,X_org(24,:))
title('original channel 24 - 20dB')
subplot(3,1,3)
plot(t,X_noisy(24,:))
title('noisy channel 24 - 20dB')