clc; clear all; close all;

load('ElecPosXYZ') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;

Radius = ModelParams.R(3) ;
load('ElecPosXYZ') ;
ElectrodePos = [] ;
Label =[];
for i=1:21
    A = ElecPos{i};
    Label{i} = num2str(A.Name);
    ElectrodePos(i,:) = Radius*A.XYZ ;
end

t = 0:0.01:2;
f = 5;
Signal = 10*sin(2*pi*f*t);

%% j
[r,c] = size(LocMat);
%index_15_dipole = [1089 1090 1091 1092 1093 1202 1203 1204 1205 1206 951 952 953 954 955];
% index_15_dipole = 37:45;
% a = 50:55;
% index_15_dipole = [index_15_dipole,a];
index_15_dipole = [204 205 206 207 208 216 217 218 219 220 229 230 231 232 233];
for i=1:15
    jahat_15_dipole(i,:) = LocMat(:,index_15_dipole(i)) ./ norm(LocMat(:,index_15_dipole(i)));
end

figure(8)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'o')
hold on
scatter3(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),'r*')
text(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),Label)

for i=1:15
    hold on
    scatter3(LocMat(1,index_15_dipole(i)),LocMat(2,index_15_dipole(i)),LocMat(3,index_15_dipole(i)),'k*')

    hold on
    plot3([LocMat(1,index_15_dipole(i)) LocMat(1,index_15_dipole(i))+jahat_15_dipole(i,1)],...
        [LocMat(2,index_15_dipole(i)) LocMat(2,index_15_dipole(i))+jahat_15_dipole(i,2)],...
        [LocMat(3,index_15_dipole(i)) LocMat(3,index_15_dipole(i))+jahat_15_dipole(i,3)],'g')
end

Q_15dipole = zeros(3951,length(t));
for i=1:15
    Q_15dipole(3*index_15_dipole(i)-2:3*index_15_dipole(i),:) = jahat_15_dipole(i,:)'*Signal;
end

M = GainMat*Q_15dipole;

for i=1:21
    figure(9)
    subplot(7,3,i)
    plot(t,M(i,:))
    title(Label{i})
end

Amp = max(M');
Amp_normalize = (Amp-mean(Amp))/std(Amp);
figure(10)
Display_Potential_3D(Radius,Amp_normalize)

%% MNE
I21 = eye(21);
alpha = 0.5;
Q_MNE = GainMat'*inv(GainMat*GainMat'+alpha*I21)*M;
A_MNE = max(Q_MNE');
for i=1:1317
    abs_MNE(i) = (sum(A_MNE(3*i-2:3*i).^2))^0.5;
end
[Amp_MNE,S_MNE] = sort(abs_MNE);
for i=1:1317
    jahat_MNE(i,:) = A_MNE(3*i-2:3*i)/norm(A_MNE(3*i-2:3*i));
end

a = (Q_15dipole-Q_MNE).^2;
RMSE_MNE = sum(a,'all')/(3951*length(t))

T = zeros(1,1317);
T(index_15_dipole) = 1;

Thr_mne = 0:0.01:floor(max(abs_MNE))+1;
S_ROC_MNE = 0;
for i=1:length(Thr_mne)
    T_mne = zeros(1,1317);
    for j=1:1317
        if (abs_MNE(j)>Thr_mne(i))
            T_mne(j) = 1;
        end
    end
    C = confusionmat(T,T_mne);
    FPR(i) = C(1,2)/(C(1,2)+C(1,1));
    TPR(i) = C(2,2)/(C(2,2)+C(2,1));
    Acc(i) = (C(2,2)+C(1,1)) / (C(1,2)+C(1,1)+C(2,2)+C(2,1));
    DD(i) = ((FPR(i)^2)+((1-TPR(i))^2))^0.5;
    if i~=1
        S_ROC_MNE = S_ROC_MNE + TPR(i)*(FPR(i-1)-FPR(i));
    end
end

figure(11)
plot(FPR,TPR)
title('ROC-MNE')
xlabel('1-Specificity')
ylabel('Sensitivity')

[MIN,I] = min(DD);
Thr_Optimum_MNE = Thr_mne(I);
dipole_selected_MNE = [];
for j=1:1317
    if (abs_MNE(j)>Thr_Optimum_MNE)
        dipole_selected_MNE = [dipole_selected_MNE,j];
    end
end
Acc__Optimum_MNE = Acc(I)
FPR__Optimum_MNE = FPR(I)
TPR__Optimum_MNE = TPR(I)
S_ROC_MNE


%% WMNE
I3 = eye(3);
omega = zeros(1317);
for i=1:1317 
    L = 0;
    for n=1:21
        L = L + GainMat(n,3*i-2:3*i)*GainMat(n,3*i-2:3*i)';
    end
    omega(i,i) = L^0.5;
end
W = zeros(3951,3951);
for i=1:3951
    z = floor((i-1)/3)+1;
    W(i,i) = omega(z,z);
end
Q_WMNE = inv(W'*W)*GainMat'*inv(GainMat*inv(W'*W)*GainMat'+alpha*I21)*M;
A_WMNE = max(Q_WMNE');
for i=1:1317
    abs_WMNE(i) = (sum(A_WMNE(3*i-2:3*i).^2))^0.5;
end
[Amp_WMNE,S_WMNE] = sort(abs_WMNE);
for i=1:1317
    jahat_WMNE(i,:) = A_WMNE(3*i-2:3*i)/norm(A_WMNE(3*i-2:3*i));
end
b = (Q_WMNE-Q_15dipole).^2;
RMSE_WMNE = sum(b,'all')/(3951*length(t));

Thr_wmne = 0:0.01:floor(max(abs_WMNE))+1;
S_ROC_WMNE = 0;
for i=1:length(Thr_wmne)
    T_wmne = zeros(1,1317);
    for j=1:1317
        if (abs_WMNE(j)>Thr_wmne(i))
            T_wmne(j) = 1;
        end
    end
    C = confusionmat(T,T_wmne);
    FPR(i) = C(1,2)/(C(1,2)+C(1,1));
    TPR(i) = C(2,2)/(C(2,2)+C(2,1));
    Acc(i) = (C(2,2)+C(1,1)) / (C(1,2)+C(1,1)+C(2,2)+C(2,1));
    DD(i) = ((FPR(i)^2)+((1-TPR(i))^2))^0.5;
    if i~=1
        S_ROC_WMNE = S_ROC_WMNE + TPR(i)*(FPR(i-1)-FPR(i));
    end
end

figure(12)
plot(FPR,TPR)
title('ROC-WMNE')
xlabel('1-Specificity')
ylabel('Sensitivity')

[MIN,I] = min(DD);
Thr_Optimum_WMNE = Thr_wmne(I);
dipole_selected_WMNE = [];
for j=1:1317
    if (abs_WMNE(j)>Thr_Optimum_WMNE)
        dipole_selected_WMNE = [dipole_selected_WMNE,j];
    end
end
Acc__Optimum_WMNE = Acc(I)
FPR__Optimum_WMNE = FPR(I)
TPR__Optimum_WMNE = TPR(I)
S_ROC_WMNE


%% LORETA
d = Resolution;
p = 1317;
A_1 = zeros(p);
for i=1:p
    for j=1:p
        L = norm(LocMat(:,i)-LocMat(:,j));
        if L==d
            A_1(i,j) = 1/6;
        end
    end
end
I1317 = eye(p);
A_0 = inv(diag(A_1 * ones(p,1))) * A_1;
AA = zeros(3*p);
for i=1:p
    for j=1:p
        AA(3*i-2:3*i,3*j-2:3*j) = A_0(i,j)*I3;
    end
end
B = (6/(d*d))*(AA-eye(3*p));


W_loreta = W * B' * B * W;
Q_loreta = inv(W_loreta'*W_loreta)*GainMat'*inv(GainMat*inv(W_loreta'*W_loreta)*GainMat'+alpha*I21)*M;

A_loreta = max(Q_loreta');
for i=1:1317
    abs_loreta(i) = (sum(A_loreta(3*i-2:3*i).^2))^0.5;
end
[Amp_loreta,S_loreta] = sort(abs_loreta);
for i=1:1317
    jahat_loreta(i,:) = A_loreta(3*i-2:3*i)/norm(A_loreta(3*i-2:3*i));
end
c = (Q_loreta-Q_15dipole).^2;
RMSE_loreta = sum(c,'all')/(3951*length(t));

Thr_loreta = 0:0.01:floor(max(abs_loreta))+1;
S_ROC_loreta = 0;
for i=1:length(Thr_loreta)
    T_loreta = zeros(1,1317);
    for j=1:1317
        if (abs_loreta(j)>Thr_loreta(i))
            T_loreta(j) = 1;
        end
    end
    C = confusionmat(T,T_loreta);
    FPR(i) = C(1,2)/(C(1,2)+C(1,1));
    TPR(i) = C(2,2)/(C(2,2)+C(2,1));
    Acc(i) = (C(2,2)+C(1,1)) / (C(1,2)+C(1,1)+C(2,2)+C(2,1));
    DD(i) = ((FPR(i)^2)+((1-TPR(i))^2))^0.5;
    if i~=1
        S_ROC_loreta = S_ROC_loreta + TPR(i)*(FPR(i-1)-FPR(i));
    end
end

figure(13)
plot(FPR,TPR)
title('ROC-loreta')
xlabel('1-Specificity')
ylabel('Sensitivity')

[MIN,I] = min(DD);
Thr_Optimum_loreta = Thr_loreta(I);
dipole_selected_loreta = [];
for j=1:1317
    if (abs_loreta(j)>Thr_Optimum_loreta)
        dipole_selected_loreta = [dipole_selected_loreta,j];
    end
end
Acc__Optimum_loreta = Acc(I)
FPR__Optimum_loreta = FPR(I)
TPR__Optimum_loreta = TPR(I)
S_ROC_loreta