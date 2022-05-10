% HW3 - 97200216
clc; clear all; close all;
%%
load('HW3_Ex1_data.mat')
fs = 100;
t = (0:10000-1)/fs;
% plot(t,X_org)
N = length(X_org);
%% 
cov_X = cov(transpose(X_org));
[U,L]=eig(cov_X);
D = (L^-0.5)*transpose(U);
X_white = D*(X_org-repmat(mean(X_org,2),1,N));

L = N/(4*fs);

w = ones(8,1);
for i=1:100
    s1 = w'*X_white;
    S = zeros(1,4*fs);
    for l=1:L
        S = S + s1(:,1+4*fs*(l-1):4*fs*l);
    end
    S = S./L;
    S_pluse = repmat(S,1,L);
    
    w_pluse = X_white*S_pluse';
    
    w = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end

s1 = w'*X_white;
%s1(1:7,:) = 0;
X1_hat = inv(D)*w*s1 + repmat(mean(X_org,2),1,N);
error_X1 = (sum(sum((X1-X1_hat).^2))/sum(sum(X1.^2)))^0.5;
for i=1:8
    figure(1)
    subplot(4,2,i)
    plot(t,X1_hat(i,:),t,X1(i,:))
    title(['channel',num2str(i)])
    legend('X1-hat','X1')
end

%%
for i=3*fs:7*fs
    cor1(i-3*fs+1) = corr2(X_org(1,1:end-i),X_org(1,i+1:end));
    cor2(i-3*fs+1) = corr2(X_org(2,1:end-i),X_org(2,i+1:end));
    cor3(i-3*fs+1) = corr2(X_org(3,1:end-i),X_org(3,i+1:end));
    cor4(i-3*fs+1) = corr2(X_org(4,1:end-i),X_org(4,i+1:end));
    cor5(i-3*fs+1) = corr2(X_org(5,1:end-i),X_org(5,i+1:end));
    cor6(i-3*fs+1) = corr2(X_org(6,1:end-i),X_org(6,i+1:end));
    cor7(i-3*fs+1) = corr2(X_org(7,1:end-i),X_org(7,i+1:end));
    cor8(i-3*fs+1) = corr2(X_org(8,1:end-i),X_org(8,i+1:end));  
end
Corr =  cor1 + cor2 + cor3 + cor4 +  cor5 + cor6 + cor7 + cor8; 
[M(1),I(1)] = max(cor1);
[M(2),I(2)] = max(cor2);
[M(3),I(3)] = max(cor3);
[M(4),I(4)] = max(cor4);
[M(5),I(5)] = max(cor5);
[M(6),I(6)] = max(cor6);
[M(7),I(7)] = max(cor7);
[M(8),I(8)] = max(cor8);
[s,i1] = sort(M);
p = (I(i1(8))+I(i1(7))+I(i1(6)))/3;
T_period = round((p+3*fs-1))/fs;
N_T = round(T_period*fs);

L = round(N/(N_T));

w = ones(8,1);
for i=1:100
    s1 = w'*X_white;
    S = zeros(1,N_T);
    for l=1:L
        S = S + s1(:,1+N_T*(l-1):N_T*l);
    end
    S = S./L;
    S_pluse = repmat(S,1,L);
    w_pluse = X_white(:,1:N_T*L)*S_pluse';
    w = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end

s1 = w'*X_white;
X1_hat_T = inv(D)*w*s1 + repmat(mean(X_org,2),1,N);
error_X1_T = (sum(sum((X1-X1_hat_T).^2))/sum(sum(X1.^2)))^0.5;
for i=1:8
    figure(2)
    subplot(4,2,i)
    plot(t,X1_hat_T(i,:),t,X1(i,:))
    title(['channel',num2str(i)])
    legend('X1-hat','X1')
end
%%
w = ones(8,1);
for i=1:100
    s1 = w'*X_white;
    S_pluse = s1.*T1;
    w_pluse = X_white*S_pluse';
    w = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end

s1 = w'*X_white;
X2_hat_T1 = inv(D)*w*s1 + repmat(mean(X_org,2),1,N);
error_X2_T1 = (sum(sum((X2-X2_hat_T1).^2))/sum(sum(X2.^2)))^0.5;
for i=1:8
    figure(3)
    subplot(4,2,i)
    plot(t,X2_hat_T1(i,:),t,X2(i,:))
    title(['channel',num2str(i)])
    legend('X2-hat','X2')
end
%%
w = ones(8,1);
for i=1:100
    s1 = w'*X_white;
    S_pluse = s1.*T2;
    w_pluse = X_white*S_pluse';
    w = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end

s1 = w'*X_white;
X2_hat_T2 = inv(D)*w*s1 + repmat(mean(X_org,2),1,N);
error_X2_T2 = (sum(sum((X2-X2_hat_T2).^2))/sum(sum(X2.^2)))^0.5;
for i=1:8
    figure(4)
    subplot(4,2,i)
    plot(t,X2_hat_T2(i,:),t,X2(i,:))
    title(['channel',num2str(i)])
    legend('X2-hat','X2')
end

%%
Z = fft(X_white');
Z = Z';
f = fs*(0:(N-1))/N;

for i=1:length(Z)
    if(f(i)>=10 && f(i)<=15)
        H(i) = 1;
    else if (f(i)>=85 && f(i)<=90)
        H(i) = 1;
        else 
            H(i) = 0;
        end
    end
end
w(:,1) = ones(8,1);
for i=1:100
    s1 = w(:,i)'*Z;
    S_pluse = s1.*H;
    w_pluse = Z*S_pluse';
    w(:,i+1) = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end
s1 = w(:,101)'*X_white;
X3_hat = inv(D)*w(:,101)*s1 + repmat(mean(X_org,2),1,N);
error_X3 = (sum(sum((X3-X3_hat).^2))/sum(sum(X3.^2)))^0.5;
Y3_hat = (fft(X3_hat'));
Y3_hat = Y3_hat';
Y3_hat = abs(Y3_hat);
for i=1:8
    figure(5)
    subplot(4,2,i)
    plot(f,Y3_hat(i,:))
    title(['channel',num2str(i)])
   % legend('X3-hat','X3')
end

%%

for i=1:100
    for j=1:8
        BP(j,i) = bandpower(X_org(j,:)-repmat(mean(X_org(j,:)),1,N),fs,[5+(i-1)*0.2 5+i*0.2]);
    end
    
end
% for i=1:8
%     figure(10)
%     subplot(4,2,i)
%     plot(BP(i,:))
% end

f_on = zeros(N,1);
for i=1:8
    a = find(abs(BP(i,:))>(10^-3));
    f_on(a) = 1;
end
ff = find(f_on==1);
f_min = 5 + ff(1)*0.2;
f_max = 5 + ff(end)*0.2;

a1 = zeros(8,8);
m1 = 0;
for i=1:N
    if(f(i)>=f_min && f(i)<=f_max)
        a1 = a1 + (Y(:,i)-mean(Y,2))*transpose(Y(:,i)-mean(Y,2));
        m1 = m1+1;
    end
    if(f(i)>=(100-f_max) && f(i)<=(100-f_min))
        a1 = a1 + (Y(:,i))*transpose(Y(:,i));
        m1 = m1+1;
    end
end
Px_FFT1 = a1 ./ m1;

[V_FFT1,D_FFT1] = eig(Px_FFT1,Cx_FFT1);
S_FFT1 = V_FFT1' * X_org;
S_denoised_FFT1 = S_FFT1;
S_denoised_FFT1(1:7,:) = 0;
X3_hat1 = inv(V_FFT1') * S_denoised_FFT1 + repmat(mean(X_org,2),1,N);
error_X3_1 = (sum(sum((X3-X3_hat1).^2))/sum(sum(X3.^2)))^0.5;

Y3_hat = (fft(X3_hat1'));
Y3_hat = Y3_hat';
Y3_hat = abs(Y3_hat);
Y3 = (fft(X3'));
Y3 = Y3';
Y3 = abs(Y3);

f = fs*(0:(N-1))/N;
for i=1:8
    figure(6)
    subplot(4,2,i)
    plot(f,Y3_hat(i,:))
    title(['channel',num2str(i)])
%     legend('X3-hat','X3')
end


%%
for i=1:length(Z)
    if(f(i)>=5 && f(i)<=25)
        H1(i) = 1;
    else if(f(i)>=75 && f(i)<=95)
        H1(i) = 1;
        else
            H1(i) = 0;
        end
    end
end
w(:,1) = ones(8,1);
for i=1:100
    s1 = w(:,i)'*Z;
    S_pluse = s1.*H1;
    w_pluse = Z*S_pluse';
    w(:,i+1) = w_pluse ./ ((sum(w_pluse.^2))^0.5);
end
s1 = w(:,101)'*X_white;
X3_hat_1 = inv(D)*w(:,101)*s1 + repmat(mean(X_org,2),1,N);
error_X3_1 = (sum(sum((X3-X3_hat_1).^2))/sum(sum(X3.^2)))^0.5;

Y3_hat = (fft(X3_hat_1'));
Y3_hat = Y3_hat';
Y3_hat = abs(Y3_hat);
for i=1:8
    figure(7)
    subplot(4,2,i)
    plot(f,Y3_hat(i,:))
    title(['channel',num2str(i)])
    %legend('X3-hat','X3')
end
