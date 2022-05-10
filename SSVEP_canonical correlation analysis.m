clc; clear all; close all;

%%
fs = 256;
t = 2:1/fs:15-1/fs;

Y_8Hz = [sin(2*pi*8*t);cos(2*pi*8*t);sin(2*pi*16*t);cos(2*pi*16*t);sin(2*pi*24*t);cos(2*pi*24*t);...
    sin(2*pi*32*t);cos(2*pi*32*t);sin(2*pi*40*t);cos(2*pi*40*t)];
Y_14Hz = [sin(2*pi*14*t);cos(2*pi*14*t);sin(2*pi*28*t);cos(2*pi*28*t);sin(2*pi*42*t);cos(2*pi*42*t)];
Y_28Hz = [sin(2*pi*28*t);cos(2*pi*28*t)];

%% 8Hz
fpass = [4 45];
load ('SSVEP_8Hz_Trial1_SUBJ2.mat')
Select = [28,42,56,63,77,84,91];
EEG = EEGSignal(Select,:);
SSVEP_8Hz_Trial1 = bandpass(EEG',fpass,fs)';
SSVEP_8Hz_Trial1(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_8Hz_Trial1',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_8Hz_Trial1',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_8Hz_Trial1',Y_28Hz');
maxR(1,:) = [r1(1),r2(1),r3(1)];
[m_R(1,1),I(1,1)] = max(maxR(1,:));
maxB(1,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(1,1),IB(1,1)] = max(maxB(1,:));

load ('SSVEP_8Hz_Trial2_SUBJ2.mat')
Select = [14,63,77,84,91];
EEG = EEGSignal(Select,:);
SSVEP_8Hz_Trial2 = bandpass(EEG',fpass,fs)';
SSVEP_8Hz_Trial2(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_8Hz_Trial2',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_8Hz_Trial2',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_8Hz_Trial2',Y_28Hz');
maxR(2,:) = [r1(1),r2(1),r3(1)];
[m_R(1,2),I(1,2)] = max(maxR(2,:));
maxB(2,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(1,2),IB(1,2)] = max(maxB(2,:));

load ('SSVEP_8Hz_Trial3_SUBJ2.mat')
Select = [24,28,30,77,90,91,105];
EEG = EEGSignal(Select,:);
SSVEP_8Hz_Trial3 = bandpass(EEG',fpass,fs)';
SSVEP_8Hz_Trial3(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_8Hz_Trial3',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_8Hz_Trial3',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_8Hz_Trial3',Y_28Hz');
maxR(3,:) = [r1(1),r2(1),r3(1)];
[m_R(1,3),I(1,3)] = max(maxR(3,:));
maxB(3,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(1,3),IB(1,3)] = max(maxB(3,:));

load ('SSVEP_8Hz_Trial4_SUBJ2.mat')
Select = [14,35,49,63,77,91];
EEG = EEGSignal(Select,:);
SSVEP_8Hz_Trial4 = bandpass(EEG',fpass,fs)';
SSVEP_8Hz_Trial4(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_8Hz_Trial4',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_8Hz_Trial4',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_8Hz_Trial4',Y_28Hz');
maxR(4,:) = [r1(1),r2(1),r3(1)];
[m_R(1,4),I(1,4)] = max(maxR(4,:));
maxB(4,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(1,4),IB(1,4)] = max(maxB(4,:));

load ('SSVEP_8Hz_Trial5_SUBJ2.mat')
Select = [24,25,63];
EEG = EEGSignal(Select,:);
SSVEP_8Hz_Trial5 = bandpass(EEG',fpass,fs)';
SSVEP_8Hz_Trial5(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_8Hz_Trial5',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_8Hz_Trial5',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_8Hz_Trial5',Y_28Hz');
maxR(5,:) = [r1(1),r2(1),r3(1)];
[m_R(1,5),I(1,5)] = max(maxR(5,:));
maxB(5,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(1,5),IB(1,5)] = max(maxB(5,:));

%% 14Hz

load ('SSVEP_14Hz_Trial1_SUBJ2.mat')
Select = [15,24,25,28];
EEG = EEGSignal(Select,:);
SSVEP_14Hz_Trial1 = bandpass(EEG',fpass,fs)';
SSVEP_14Hz_Trial1(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_14Hz_Trial1',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_14Hz_Trial1',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_14Hz_Trial1',Y_28Hz');
maxR(6,:) = [r1(1),r2(1),r3(1)];
[m_R(2,1),I(2,1)] = max(maxR(6,:));
maxB(6,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(2,1),IB(2,1)] = max(maxB(6,:));

load ('SSVEP_14Hz_Trial2_SUBJ2.mat')
Select = [14,15,16,21,23,24,25];
EEG = EEGSignal(Select,:);
SSVEP_14Hz_Trial2 = bandpass(EEG',fpass,fs)';
SSVEP_14Hz_Trial2(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_14Hz_Trial2',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_14Hz_Trial2',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_14Hz_Trial2',Y_28Hz');
maxR(7,:) = [r1(1),r2(1),r3(1)];
[m_R(2,2),I(2,2)] = max(maxR(7,:));
maxB(7,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(2,2),IB(2,2)] = max(maxB(7,:));

load ('SSVEP_14Hz_Trial3_SUBJ2.mat')
Select = [24,25,28,30];
EEG = EEGSignal(Select,:);
SSVEP_14Hz_Trial3 = bandpass(EEG',fpass,fs)';
SSVEP_14Hz_Trial3(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_14Hz_Trial3',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_14Hz_Trial3',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_14Hz_Trial3',Y_28Hz');
maxR(8,:) = [r1(1),r2(1),r3(1)];
[m_R(2,3),I(2,3)] = max(maxR(8,:));
maxB(8,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(2,3),IB(2,3)] = max(maxB(8,:));

load ('SSVEP_14Hz_Trial4_SUBJ2.mat')
Select = [6,24,26,27,28,53];
EEG = EEGSignal(Select,:);
SSVEP_14Hz_Trial4 = bandpass(EEG',fpass,fs)';
SSVEP_14Hz_Trial4(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_14Hz_Trial4',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_14Hz_Trial4',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_14Hz_Trial4',Y_28Hz');
maxR(9,:) = [r1(1),r2(1),r3(1)];
[m_R(2,4),I(2,4)] = max(maxR(9,:));
maxB(9,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(2,4),IB(2,4)] = max(maxB(9,:));

load ('SSVEP_14Hz_Trial5_SUBJ2.mat')
Select = [18,46,49,55];
EEG = EEGSignal(Select,:);
SSVEP_14Hz_Trial5 = bandpass(EEG',fpass,fs)';
SSVEP_14Hz_Trial5(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_14Hz_Trial5',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_14Hz_Trial5',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_14Hz_Trial5',Y_28Hz');
maxR(10,:) = [r1(1),r2(1),r3(1)];
[m_R(2,5),I(2,5)] = max(maxR(10,:));
maxB(10,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(2,5),IB(2,5)] = max(maxB(10,:));

%% 28Hz
load ('SSVEP_28Hz_Trial1_SUBJ2.mat')
Select = [18,46,49,55];
EEG = EEGSignal(Select,:);
SSVEP_28Hz_Trial1 = bandpass(EEG',fpass,fs)';
SSVEP_28Hz_Trial1(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_28Hz_Trial1',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_28Hz_Trial1',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_28Hz_Trial1',Y_28Hz');
maxR(11,:) = [r1(1),r2(1),r3(1)];
[m_R(3,1),I(3,1)] = max(maxR(11,:));
maxB(11,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(3,1),IB(3,1)] = max(maxB(11,:));

load ('SSVEP_28Hz_Trial2_SUBJ2.mat')
Select = [1:128];
EEG = EEGSignal(Select,:);
SSVEP_28Hz_Trial2 = bandpass(EEG',fpass,fs)';
SSVEP_28Hz_Trial2(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_28Hz_Trial2',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_28Hz_Trial2',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_28Hz_Trial2',Y_28Hz');
maxR(12,:) = [r1(1),r2(1),r3(1)];
[m_R(3,2),I(3,2)] = max(maxR(12,:));
maxB(12,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(3,2),IB(3,2)] = max(maxB(12,:));

load ('SSVEP_28Hz_Trial3_SUBJ2.mat')
Select = [18,52,112,113];
EEG = EEGSignal(Select,:);
SSVEP_28Hz_Trial3 = bandpass(EEG',fpass,fs)';
SSVEP_28Hz_Trial3(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_28Hz_Trial3',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_28Hz_Trial3',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_28Hz_Trial3',Y_28Hz');
maxR(13,:) = [r1(1),r2(1),r3(1)];
[m_R(3,3),I(3,3)] = max(maxR(13,:));
maxB(13,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(3,3),IB(3,3)] = max(maxB(13,:));

load ('SSVEP_28Hz_Trial4_SUBJ2.mat')
Select = [18,52,112,113];
EEG = EEGSignal(Select,:);
SSVEP_28Hz_Trial4 = bandpass(EEG',fpass,fs)';
SSVEP_28Hz_Trial4(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_28Hz_Trial4',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_28Hz_Trial4',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_28Hz_Trial4',Y_28Hz');
maxR(14,:) = [r1(1),r2(1),r3(1)];
[m_R(3,4),I(3,4)] = max(maxR(14,:));
maxB(14,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(3,4),IB(3,4)] = max(maxB(14,:));

load ('SSVEP_28Hz_Trial5_SUBJ2.mat')
Select = [6,15,18,20,25,50,98];
EEG = EEGSignal(Select,:);
SSVEP_28Hz_Trial5 = bandpass(EEG',fpass,fs)';
SSVEP_28Hz_Trial5(:,1:2*fs)=[];
[A1,B1,r1] = canoncorr(SSVEP_28Hz_Trial5',Y_8Hz');
[A2,B2,r2] = canoncorr(SSVEP_28Hz_Trial5',Y_14Hz');
[A3,B3,r3] = canoncorr(SSVEP_28Hz_Trial5',Y_28Hz');
maxR(15,:) = [r1(1),r2(1),r3(1)];
[m_R(3,5),I(3,5)] = max(maxR(15,:));
maxB(15,:) = [max(B1(:,1)),max(B2(:,1)),max(B3(:,1))];
[m_B(3,5),IB(3,5)] = max(maxB(15,:));

%%
for i=1:3
    for j=1:5
        if (I(i,j)==1)
            Label_with_r(i,j) = 8;
        else if (I(i,j)==2)
                Label_with_r(i,j) = 14;
            else
                Label_with_r(i,j) = 28;
            end
        end
    end
end

for i=1:3
    for j=1:5
        if (IB(i,j)==1)
            Label_with_B(i,j) = 8;
        else if (IB(i,j)==2)
                Label_with_B(i,j) = 14;
            else
                Label_with_B(i,j) = 28;
            end
        end
    end
end