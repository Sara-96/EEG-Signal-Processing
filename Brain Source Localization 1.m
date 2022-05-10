clc; clear all; close all;

%% a

load('ElecPosXYZ') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;

figure(1)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'o')

%% b
Radius = ModelParams.R(3) ;
load('ElecPosXYZ') ;
ElectrodePos = [] ;
Label =[];
for i=1:21
    A = ElecPos{i};
    Label{i} = num2str(A.Name);
    ElectrodePos(i,:) = Radius*A.XYZ ;
end

figure(2)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'o')
hold on
scatter3(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),'r*')
text(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),Label)

%% c
[r,c] = size(LocMat);
S_dipole = randi([1 c],[1 1]);
figure(3)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'o')
hold on
scatter3(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),'r*')
text(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),Label)
hold on
scatter3(LocMat(1,S_dipole),LocMat(2,S_dipole),LocMat(3,S_dipole),'k*')

jahat_dipole = LocMat(:,S_dipole) / norm(LocMat(:,S_dipole));
hold on
plot3([LocMat(1,S_dipole) LocMat(1,S_dipole)+jahat_dipole(1)],[LocMat(2,S_dipole) LocMat(2,S_dipole)+jahat_dipole(2)],...
    [LocMat(3,S_dipole) LocMat(3,S_dipole)+jahat_dipole(3)],'g')

%% d
t = 0:0.01:2;
f = 5;
Signal = 10*sin(2*pi*f*t);

G = GainMat(:,3*S_dipole-2:3*S_dipole);
Q = jahat_dipole*Signal;

M = G*Q;
  
for i=1:21
    figure(4)
    subplot(7,3,i)
    plot(t,M(i,:))
    title(Label{i})
end

%% e
Amp = max(M');
Amp_normalize = (Amp-mean(Amp))/std(Amp);
figure(5)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'ro')
hold on
SI = repmat(50,21,1);
scatter3(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),SI,Amp_normalize,'*')
colorbar
text(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),Label)
hold on
scatter3(LocMat(1,S_dipole),LocMat(2,S_dipole),LocMat(3,S_dipole),'k*')

jahat_dipole = LocMat(:,S_dipole) / norm(LocMat(:,S_dipole));
hold on
plot3([LocMat(1,S_dipole) LocMat(1,S_dipole)+jahat_dipole(1)],[LocMat(2,S_dipole) LocMat(2,S_dipole)+jahat_dipole(2)],...
    [LocMat(3,S_dipole) LocMat(3,S_dipole)+jahat_dipole(3)],'g')


figure(6)
SI = repmat(50,21,1);
scatter3(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),SI,Amp_normalize,'filled')
colorbar
text(ElectrodePos(:,1),ElectrodePos(:,2),ElectrodePos(:,3),Label)
%% f
figure(7)
Display_Potential_3D(Radius,Amp_normalize)

%% g
% MNE
Q_real = zeros(3951,length(t));
Q_real(3*S_dipole-2:3*S_dipole,:) = jahat_dipole*Signal;;

I21 = eye(21);
alpha = 0.5;
Q_MNE = GainMat'*inv(GainMat*GainMat'+alpha*I21)*M;

% WMNE
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

% LORETA
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
W_loreta = W*B'*B*W;
Q_loreta = inv(W_loreta'*W_loreta)*GainMat'*inv(GainMat*inv(W_loreta'*W_loreta)*GainMat'+alpha*I21)*M;


%% h

A_MNE = max(Q_MNE');
for i=1:1317
    Amp_MNE(i) = (sum(A_MNE(3*i-2:3*i).^2))^0.5;
end
[Amp_dipole_MNE,S_dipole_MNE] = max(Amp_MNE)
jahat_MNE = A_MNE(3*S_dipole_MNE-2:3*S_dipole_MNE)'/Amp_dipole_MNE

A_WMNE = max(Q_WMNE');
for i=1:1317
    Amp_WMNE(i) = (sum(A_WMNE(3*i-2:3*i).^2))^0.5;
end
[Amp_dipole_WMNE,S_dipole_WMNE] = max(Amp_WMNE)
jahat_WMNE = A_WMNE(3*S_dipole_WMNE-2:3*S_dipole_WMNE)'/Amp_dipole_WMNE

A_loreta = max(Q_loreta');
for i=1:1317
    Amp_loreta(i) = (sum(A_loreta(3*i-2:3*i).^2))^0.5;
end
[Amp_dipole_loreta,S_dipole_loreta] = max(Amp_loreta)
jahat_loreta = A_loreta(3*S_dipole_loreta-2:3*S_dipole_loreta)'/Amp_dipole_loreta

%% i
a = (Q_real-Q_MNE).^2;
RMSE_MNE = sum(a,'all')/(3951*length(t))
d_MNE = norm(LocMat(:,S_dipole_MNE) - LocMat(:,S_dipole))
Djahat_MNE = jahat_MNE - jahat_dipole

b = (Q_real-Q_WMNE).^2;
RMSE_WMNE = sum(b,'all')/(3951*length(t))
d_WMNE = norm(LocMat(:,S_dipole_WMNE) - LocMat(:,S_dipole))
Djahat_WMNE = jahat_WMNE - jahat_dipole

c = (Q_real-Q_loreta).^2;
RMSE_loreta = sum(c,'all')/(3951*length(t))
d_loreta = norm(LocMat(:,S_dipole_loreta) - LocMat(:,S_dipole))
Djahat_loreta = jahat_loreta - jahat_dipole

