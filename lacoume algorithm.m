clc
clear all;
close all;
load('Ex2.mat')
%%
A = rand(2);
X = A*S;
[r,N] = size(X);
n = 1:N;
[S_prime] = Lacoume(X,A,200);

CORR(1,1) = corr2(S_prime(1,:),S(1,:));
CORR(2,1) = corr2(S_prime(2,:),S(1,:));
CORR(1,2) = corr2(S_prime(1,:),S(2,:));
CORR(2,2) = corr2(S_prime(2,:),S(2,:));

[m1,i1] = max(abs(CORR(:,1)));
[m2,i2] = max(abs(CORR(:,2)));

if (i1~=i2)
    y11 = S_prime(i1,:);
    y22 = S_prime(i2,:);
else
    [m3,i3] = max(abs(CORR(1,:)));
    [m4,i4] = max(abs(CORR(2,:)));
    y11 = S_prime(i3,:);
    y22 = S_prime(i4,:);
end

p1 = mean(S(1,:)) / mean(y11);
p2 = mean(S(2,:)) / mean(y22);

y1 = p1*y11;
y2 = p2*y22;
e1 = sum((y1-S(1,:)).^2)/N;
e2 = sum((y2-S(2,:)).^2)/N;

correlation = [m1 m2]
square_error = [e1 e2]

figure(1)
plot(n,S(1,:),n,y1)
title('channel 1')
legend('real','after seperation')

figure(2)
plot(n,S(2,:),n,y2)
title('channel 2')
legend('real','after seperation')
%%

for i = 1:5
    Var(i) = 0.1*i;
    noise = random('Normal',0,Var(i),[r N]);
    XwithNoise = X + noise;
    [S_noise] = Lacoume(XwithNoise,A,40-5*i);

    CORRnoise(1,1) = corr2(S_noise(1,:),S(1,:));
    CORRnoise(2,1) = corr2(S_noise(2,:),S(1,:));
    CORRnoise(1,2) = corr2(S_noise(1,:),S(2,:));
    CORRnoise(2,2) = corr2(S_noise(2,:),S(2,:));

    [m1,i1] = max(abs(CORRnoise(:,1)));
    [m2,i2] = max(abs(CORRnoise(:,2)));
    
    if (i1~=i2)
        y11 = S_noise(i1,:);
        y22 = S_noise(i2,:);
        C1(i) = m1;
        C2(i) = m2;
    else
        [m3,i3] = max(abs(CORRnoise(1,:)));
        [m4,i4] = max(abs(CORRnoise(2,:)));
        y11 = S_noise(i3,:);
        y22 = S_noise(i4,:);
        C1(i) = m3;
        C2(i) = m4;
    end

    p1 = mean(S(1,:)) / mean(y11);
    p2 = mean(S(2,:)) / mean(y22);

    y1 = p1*y11;
    y2 = p2*y22;
    e1noise(i) = sum((y1-S(1,:)).^2)/N;
    e2noise(i) = sum((y2-S(2,:)).^2)/N;
    ss = sum(X(1,:).^2) + sum(X(2,:).^2);
    nn = sum(noise(1,:).^2) + sum(noise(2,:).^2);
    SNR(i) = 10*log10(ss./nn);
end

figure(3)
subplot(2,1,1)
plot(SNR,C1)
title('channel 1')
xlabel('SNR')
ylabel('correlation')
subplot(2,1,2)
plot(SNR,C2)
title('channel 2')
xlabel('SNR')
ylabel('correlation')

figure(4)
subplot(2,1,1)
plot(SNR,e1noise)
title('channel 1')
xlabel('SNR')
ylabel('square error')
subplot(2,1,2)
plot(SNR,e2noise)
title('channel 2')
xlabel('SNR')
ylabel('square error')
    