clc;
close all; clear all;

load ('ERP_EEG.mat')
t = (0:239)/240;
N = 100:100:2500;
sum_ERP = zeros(240,1);
for i=1:length(N)
    sum_ERP = sum_ERP + ERP_EEG(:,N(i));
    M(:,i) = sum_ERP ./ i;
    figure(1)
    subplot(5,5,i)
    plot(t,M(:,i))
end

for i=1:length(N)
    MAX(i) = max(abs(M(:,i)));
    n(i) = i;
end
figure(2)
title('max')
plot(n,MAX)

for i=2:length(N)
    error(i) = rms(M(:,i)-M(:,i-1));
end
figure(3)
title('root mean square')
plot(error)

%%
sum_ERP1 = zeros(240,1);
for i=1:200
    sum_ERP1 = sum_ERP1 + ERP_EEG(:,i);
end
mean200 = sum_ERP1 / 200;

LAG = 100*ones(1,200);
mean200_new = mean200;
for m=1:1000
%     max(abs(LAG(i)))
    sum_ERP2 = zeros(240,1);
    for i=1:200
        [acor,lag] = xcorr(mean200_new,ERP_EEG(:,i));
        [~,I] = max(acor);
        LAG(i) = lag(I);
        if (LAG(i)>0)
            ERP_NEW = [zeros(abs(LAG(i)),1);ERP_EEG(1:end-abs(LAG(i)),i)];
        else if (LAG(i)<0)
                ERP_NEW = [ERP_EEG(abs(LAG(i))+1:end,i);zeros(abs(LAG(i)),1)];
            else if (LAG(i)==0)
                ERP_NEW = ERP_EEG(:,i);
            end
            end
        end
        sum_ERP2 = sum_ERP2 + ERP_NEW;
    end
    mean200_new = sum_ERP2/200;
end    


sum_ERP3 = zeros(240,1);
for i=1:2550
    sum_ERP3 = sum_ERP3 + ERP_EEG(:,i);
end
mean2550 = sum_ERP3 / 2550;


figure(4)
plot(t,[mean200,mean200_new,mean2550])
xlabel('time(s)')
legend('Mean200','Mean200 with lag','Mean 2550')