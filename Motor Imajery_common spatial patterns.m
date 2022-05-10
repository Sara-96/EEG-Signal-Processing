clc; clear all; close all;

%%
load('AllData_Feet_Hand.mat')
load('AllElectrodes.mat')

fs = 256;
class0 = find(TrainLabel == 0);
class1 = find(TrainLabel == 1);

Mean = mean(TrainData,2);
TrainData_normalized = TrainData - repmat(Mean,1,256,1);

C0 = zeros(30,30);
for m=1:length(class0)
    a = TrainData_normalized(:,:,class0(m));
    C0 = C0 + a*a'/256;
end
C0 = C0/length(class0);

C1 = zeros(30,30);
for m=1:length(class1)
    a = TrainData_normalized(:,:,class1(m));
    C1 = C1 + a*a'/256;
end
C1 = C1/length(class1);

[V,D,W] = eig(C0,C1);
F = 1;
W_csp = [W(:,1:F),W(:,end-F+1:end)];

Y_csp_train_kol0 = [];
Y_csp_train_kol1 = [];
for i=1:180
    Y_csp_train(:,:,i) = W_csp'*TrainData_normalized(:,:,i);
    if TrainLabel(i)==0
        Y_csp_train_kol0 = [Y_csp_train_kol0,Y_csp_train(:,:,i)];
    else
        Y_csp_train_kol1 = [Y_csp_train_kol1,Y_csp_train(:,:,i)];
    end
end


figure(1)
subplot(2,1,1)
plot(Y_csp_train_kol0(1,:))
title('Bottom CSP filter - hand')
subplot(2,1,2)
plot(Y_csp_train_kol0(2,:))
title('Top CSP filter - hand')
    
figure(2)
subplot(2,1,1)
plot(Y_csp_train_kol1(1,:))
title('Top CSP filter - feet')
subplot(2,1,2)
plot(Y_csp_train_kol1(2,:))
title('Bottom CSP filter - feet')


SS0(:,:) = var(Y_csp_train(:,:,class1),1,2);
SS1(:,:) = var(Y_csp_train(:,:,class0),1,2);

Mean2 = mean(TestData,2);
TestData_normalized = TestData - repmat(Mean2,1,256,1);
for i=1:60
    Y_csp_test(:,:,i) = W_csp'*TestData_normalized(:,:,i);
end

%%
[elocsX,elocsY] = pol2cart(pi/180*[AllElectrodes.theta],[AllElectrodes.radius]);
E = [37,3,7,38,40,42,10,45,47,15,48,50,13,52,18,55,32,20,21,22,23,57,58,59,60,31,26,27,63,64];
elocsX1 = elocsX(E);
elocsY1 = elocsY(E);
LL = char(AllElectrodes.labels);
figure(3)
plottopomap(elocsX1', elocsY1', LL(E)', W_csp(:,1));
figure(4)
plottopomap(elocsX1', elocsY1', LL(E)', W_csp(:,2));

%%
for i=1:4
    Test = TrainData(:,:,45*(i-1)+1:45*i);
    Train = TrainData;
    Train(:,:,45*(i-1)+1:45*i) = [];
    TLabel = TrainLabel;
    TLabel(45*(i-1)+1:45*i) = [];
    cl0 = find(TLabel ==0);
    cl1 = find(TLabel ==1);

    Mean3 = mean(Train,2);
    Train_N = Train - repmat(Mean3,1,256,1);

    C0 = zeros(30,30);
    for m=1:length(cl0)
        a = TrainData_normalized(:,:,cl0(m));
        C0 = C0 + a*a'/256;
    end
    C0 = C0/length(cl0);

    C1 = zeros(30,30);
    for m=1:length(cl1)
        a = TrainData_normalized(:,:,cl1(m));
        C1 = C1 + a*a'/256;
    end
    C1 = C1/length(cl1);

    [V,D,W1] = eig(C0,C1);
    for F=1:15
        Y_csp_Train_4fold=[]; 
        Features_train=[];
        Y_csp_Test_4fold=[];
        Features_test=[];
        
        W_csp_4fold = [W1(:,1:F),W1(:,end-F+1:end)];
        for k=1:135
            Y_csp_Train_4fold(:,:,k) = W_csp_4fold'*Train_N(:,:,k);
        end
        Features_train(:,:) = var(Y_csp_Train_4fold(:,:,:),1,2);
        
        Mean4 = mean(Test,2);
        Test_N = Test - repmat(Mean4,1,256,1);
        for k=1:45
            Y_csp_Test_4fold(:,:,k) = W_csp_4fold'*Test_N(:,:,k);
        end
        Features_test(:,:) = var(Y_csp_Test_4fold(:,:,:),1,2);
        
        svmStruct = fitcsvm(Features_train',TLabel);
        TestLabel_4fold = predict(svmStruct,transpose(Features_test)); 
        
        C = confusionmat(TrainLabel(45*(i-1)+1:45*i)',TestLabel_4fold);
        acc(i,F) = trace(C)/45*100;
         
    end
end
ACC_mean = mean(acc);
[~,F_optimom] = max(ACC_mean);

%%
Y_csp_Train_optimom=[];
Y_csp_Test_optimom=[];
Features_train_optimom=[];
Features_test_optimom=[];

W_csp_optimom = [W(:,1:F_optimom),W(:,end-F_optimom+1:end)];
for k=1:180
	Y_csp_Train_optimom(:,:,k) = W_csp_optimom'*TrainData_normalized(:,:,k);
end
Features_train_optimom(:,:) = var(Y_csp_Train_optimom,1,2);
        
for k=1:60
    Y_csp_Test_optimom(:,:,k) = W_csp_optimom'*TestData_normalized(:,:,k);
end
Features_test_optimom(:,:) = var(Y_csp_Test_optimom,1,2);

svmStruct1 = fitcsvm(transpose(Features_train_optimom),TrainLabel);
TestLabel_predict = predict(svmStruct1,transpose(Features_test_optimom));

C_test = confusionmat(TestLabel,TestLabel_predict);
ACC_test = trace(C_test)/60*100;

