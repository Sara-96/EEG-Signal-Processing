clc;
clear all;
close all;
load('Ex3.mat')

%%
cov_X = cov(transpose(X));
[U,L]=eig(cov_X);
D = (L^-0.5)*transpose(U);
Y = D*X;
cov_Y = cov(transpose(Y))

figure(1)
p = scatter3(X(1,:),X(2,:),X(3,:),'yellow')
hold on
p = plot3([0,U(1,1)],[0,U(2,1)],[0,U(3,1)])
c = p.Color;
p.Color = 'black'
hold on
p = plot3([0,U(1,2)],[0,U(2,2)],[0,U(3,2)])
c = p.Color;
p.Color = 'blue'
hold on
p = plot3([0,U(1,3)],[0,U(2,3)],[0,U(3,3)])
c = p.Color;
p.Color = 'red'
figure(2)
scatter3(Y(1,:),Y(2,:),Y(3,:))

%%
[coeff,score,latent] = pca(transpose(X));
L_prime = latent.^-0.5
D_prime = [L_prime(1).*coeff(:,1)' ;L_prime(2).*coeff(:,2)'; L_prime(3).*coeff(:,3)'];
Y_prime = D_prime*X;
cov_Y_prime = cov(transpose(Y_prime))
figure(3)
scatter3(Y_prime(1,:),Y_prime(2,:),Y_prime(3,:))
%%
[U_svd,S_svd,V_svd] = svd(cov_X);
D_svd = (S_svd^-0.5)*transpose(U_svd);
Y_svd = D_svd*X;
cov_Y_svd = cov(transpose(Y_svd))
figure(4)
scatter3(Y_svd(1,:),Y_svd(2,:),Y_svd(3,:))
