clc;
close all; clear all;

fs = 256;
load ('SSVEP_8Hz_Trial1_SUBJ2.mat')
SSVEP_8Hz_Trial1 = EEGSignal;
load ('SSVEP_8Hz_Trial2_SUBJ2.mat')
SSVEP_8Hz_Trial2 = EEGSignal;
load ('SSVEP_8Hz_Trial3_SUBJ2.mat')
SSVEP_8Hz_Trial3 = EEGSignal;
load ('SSVEP_8Hz_Trial4_SUBJ2.mat')
SSVEP_8Hz_Trial4 = EEGSignal;
load ('SSVEP_8Hz_Trial5_SUBJ2.mat')
SSVEP_8Hz_Trial5 = EEGSignal;

load ('SSVEP_14Hz_Trial1_SUBJ2.mat')
SSVEP_14Hz_Trial1 = EEGSignal;
load ('SSVEP_14Hz_Trial2_SUBJ2.mat')
SSVEP_14Hz_Trial2 = EEGSignal;
load ('SSVEP_14Hz_Trial3_SUBJ2.mat')
SSVEP_14Hz_Trial3 = EEGSignal;
load ('SSVEP_14Hz_Trial4_SUBJ2.mat')
SSVEP_14Hz_Trial4 = EEGSignal;
load ('SSVEP_14Hz_Trial5_SUBJ2.mat')
SSVEP_14Hz_Trial5 = EEGSignal;

load ('SSVEP_28Hz_Trial1_SUBJ2.mat')
SSVEP_28Hz_Trial1 = EEGSignal;
load ('SSVEP_28Hz_Trial2_SUBJ2.mat')
SSVEP_28Hz_Trial2 = EEGSignal;
load ('SSVEP_28Hz_Trial3_SUBJ2.mat')
SSVEP_28Hz_Trial3 = EEGSignal;
load ('SSVEP_28Hz_Trial4_SUBJ2.mat')
SSVEP_28Hz_Trial4 = EEGSignal;
load ('SSVEP_28Hz_Trial5_SUBJ2.mat')
SSVEP_28Hz_Trial5 = EEGSignal;

%%
t = (513:3840)/256;
s = randi([1 128],[15 9]);
%SSVEP_8Hz_Trial1
for i=1:9
    S_filtered = bandpass(SSVEP_8Hz_Trial1(s(1,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(1)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(2)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_8Hz_Trial2
for i=1:9
    S_filtered = bandpass(SSVEP_8Hz_Trial2(s(2,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(3)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(4)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_8Hz_Trial3
for i=1:9
    S_filtered = bandpass(SSVEP_8Hz_Trial3(s(3,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(5)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(6)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_8Hz_Trial4
for i=1:9
    S_filtered = bandpass(SSVEP_8Hz_Trial4(s(4,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(7)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(8)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_8Hz_Trial5
for i=1:9
    S_filtered = bandpass(SSVEP_8Hz_Trial5(s(5,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(9)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(10)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_14Hz_Trial1
for i=1:9
    S_filtered = bandpass(SSVEP_14Hz_Trial1(s(6,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(11)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(12)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_14Hz_Trial2
for i=1:9
    S_filtered = bandpass(SSVEP_14Hz_Trial2(s(7,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(13)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(14)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_14Hz_Trial3
for i=1:9
    S_filtered = bandpass(SSVEP_14Hz_Trial3(s(8,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(15)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(16)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_14Hz_Trial4
for i=1:9
    S_filtered = bandpass(SSVEP_14Hz_Trial4(s(9,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(17)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(18)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end
    
 
%SSVEP_14Hz_Trial5
for i=1:9
    S_filtered = bandpass(SSVEP_14Hz_Trial5(s(10,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(19)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(20)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_28Hz_Trial1
for i=1:9
    S_filtered = bandpass(SSVEP_28Hz_Trial1(s(11,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(21)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(22)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_28Hz_Trial2
for i=1:9
    S_filtered = bandpass(SSVEP_28Hz_Trial2(s(12,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(23)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(24)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_28Hz_Trial3
for i=1:9
    S_filtered = bandpass(SSVEP_28Hz_Trial3(s(13,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(25)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(26)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_28Hz_Trial4
for i=1:9
    S_filtered = bandpass(SSVEP_28Hz_Trial4(s(14,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(27)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(28)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end

%SSVEP_28Hz_Trial5
for i=1:9
    S_filtered = bandpass(SSVEP_28Hz_Trial5(s(15,i),:),[4 60],fs);
    S_deleted = S_filtered(513:end);
    
    figure(29)
    subplot(3,3,i)
    plot(t,S_deleted)
    title(['channel ',num2str(s(1,i))])
    xlabel('time(s)')
    
    [pxx,f] = pwelch(S_deleted,[],[],[],fs);
    figure(30)
    subplot(3,3,i)
    plot(f,pxx)
    title(['channel ',num2str(s(1,i))])
    xlabel('f(Hz)')
end
    