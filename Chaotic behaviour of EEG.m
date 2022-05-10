clc; clear all; close all;

%% logistic
N = 1:1000;
% 0<A<1
for i=1:3
    A = random('Uniform',0,1,[1 1]);
    x(1) = random('Uniform',0,1,[1 1]);
    
    for n=1:999
        x(n+1) = A*x(n)*(1-x(n));
    end
    figure(1)
    subplot(3,1,i)
    plot(N(901:1000),x(901:1000))
    title(['A=',num2str(A),' and x0=',num2str(x(1))])
end

% 1<A<3
for i=1:3
    A = random('Uniform',1,3,[1 1]);
    for j=1:2
        x(1) = random('Uniform',0,1,[1 1]);
    
        for n=1:999
            x(n+1) = A*x(n)*(1-x(n));
        end
        figure(2)
        subplot(3,2,j+2*(i-1))
        plot(N(901:1000),x(901:1000))
        title(['A=',num2str(A),' and x0=',num2str(x(1)),' and T>8'])
    end
end

% 3<A<3.449490 T=2
for i=1:3
    A = random('Uniform',3,3.449490,[1 1]);
    for j=1:2
        x(1) = random('Uniform',0,1,[1 1]);
    
        for n=1:999
            x(n+1) = A*x(n)*(1-x(n));
        end
        figure(3)
        subplot(3,2,j+2*(i-1))
        plot(N(901:1000),x(901:1000))
        title(['A=',num2str(A),' and x0=',num2str(x(1)),' and T=2'])
    end
end

% 3.449490<A<3.544090  T=4
for i=1:3
    A = random('Uniform',3.449490,3.544090,[1 1]);
    for j=1:2
        x(1) = random('Uniform',0,1,[1 1]);
    
        for n=1:999
            x(n+1) = A*x(n)*(1-x(n));
        end
        figure(4)
        subplot(3,2,j+2*(i-1))
        plot(N(901:1000),x(901:1000))
        title(['A=',num2str(A),' and x0=',num2str(x(1)),' and T=4'])
    end
end

% 3.544090<A<3.564407  T=8
for i=1:3
    A = random('Uniform',3.544090,3.564407,[1 1]);
    for j=1:2
        x(1) = random('Uniform',0,1,[1 1]);
    
        for n=1:999
            x(n+1) = A*x(n)*(1-x(n));
        end
        figure(5)
        subplot(3,2,j+2*(i-1))
        plot(N(901:1000),x(901:1000))
        title(['A=',num2str(A),' and x0=',num2str(x(1)),' and T=8'])
    end
end


% 3.564407<A<4  T>8
for i=1:3
    A = random('Uniform',3.564407,4,[1 1]);
    for j=1:3
        x(1) = random('Uniform',0,1,[1 1]);
    
        for n=1:999
            x(n+1) = A*x(n)*(1-x(n));
        end
        figure(6)
        subplot(3,3,j+3*(i-1))
        plot(N(901:1000),x(901:1000))
        title(['A=',num2str(A),' and x0=',num2str(x(1)),' and T>8'])
    end
end

%%
x = [];
P = [];
L = 3e3;
A = random('Uniform',0,4,[1 L]);
for i=1:L
    x(1,1) = random('Uniform',0,1,[1 1]);
    for n=1:999
        x(n+1,1) = A(i)*x(n,1)*(1-x(n,1));
    end
    x_inf = [repmat(A(i),20,1),x(981:1000)];
    P = [P; x_inf];
end
figure(7)
plot(P(:,1),P(:,2),'.')
title('bifurcation')
xlabel('A')
ylabel('x_inf')

P2 = [];
L = 1e4;
A = random('Uniform',2.5,4,[1 L]);
for i=1:L
    x(1,1) = random('Uniform',0,1,[1 1]);
    for n=1:999
        x(n+1,1) = A(i)*x(n,1)*(1-x(n,1));
    end
    x_inf = [repmat(A(i),10,1),x(991:1000)];
    P2 = [P2; x_inf];
end
figure(8)
plot(P2(:,1),P2(:,2),'k.')
title('bifurcation - large scale')
xlabel('A')
ylabel('x-inf')

%%
clc; clear all; close all;

%% Lorenz
b = 2.67;
s = 10;
dt = 0.01;
t = 0:dt:30;
%% r<24.74
for n=1:3
    r = random('Uniform',0,24.74,[1 1]);
    for j=1:2
        x(1) = random('Uniform',1,10,[1 1]);
        y(1) = random('Uniform',1,10,[1 1]);
        z(1) = random('Uniform',1,10,[1 1]);
        for i=1:length(t)-1
            dx = s*(y(i)-x(i))*dt;
            dy = (r*x(i)-y(i)-x(i)*z(i))*dt;
            dz = (x(i)*y(i)-b*z(i))*dt;
            
            x(i+1) = x(i) + dx;
            y(i+1) = y(i) + dy;
            z(i+1) = z(i) + dz;
        end
        figure(3*n-2)
        subplot(3,2,j)
        plot(t,x)
        title(['r=',num2str(r),',x0=',num2str(x(1)),',y0=',num2str(y(1)),',z0=',num2str(z(1))])
        xlabel('t')
        ylabel('x')
        subplot(3,2,2+j)
        plot(t,y)
        xlabel('t')
        ylabel('y')
        subplot(3,2,4+j)
        plot(t,z)
        xlabel('t')
        ylabel('z')
        
        figure(3*n-1)
        subplot(1,3,1)
        plot(x,y)
        xlabel('x')
        ylabel('y')
        subplot(1,3,2)
        plot(x,z)
        xlabel('x')
        ylabel('z')
        subplot(1,3,3)
        plot(y,z)
        xlabel('y')
        ylabel('z')
        
        figure(3*n)
        plot3(x,y,z)
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
end

%% 24.74<r<28
for n=1:3
    r = random('Uniform',24.74,28,[1 1]);
    for j=1:2
        x(1) = random('Uniform',1,10,[1 1]);
        y(1) = random('Uniform',1,10,[1 1]);
        z(1) = random('Uniform',1,10,[1 1]);
        for i=1:length(t)-1
            dx = s*(y(i)-x(i))*dt;
            dy = (r*x(i)-y(i)-x(i)*z(i))*dt;
            dz = (x(i)*y(i)-b*z(i))*dt;
            
            x(i+1) = x(i) + dx;
            y(i+1) = y(i) + dy;
            z(i+1) = z(i) + dz;
        end
        figure(3*n+7)
        subplot(3,2,j)
        plot(t,x)
        title(['r=',num2str(r),',x0=',num2str(x(1)),',y0=',num2str(y(1)),',z0=',num2str(z(1))])
        xlabel('t')
        ylabel('x')
        subplot(3,2,2+j)
        plot(t,y)
        xlabel('t')
        ylabel('y')
        subplot(3,2,4+j)
        plot(t,z)
        xlabel('t')
        ylabel('z')
        
        figure(3*n+8)
        subplot(1,3,1)
        plot(x,y)
        xlabel('x')
        ylabel('y')
        subplot(1,3,2)
        plot(x,z)
        xlabel('x')
        ylabel('z')
        subplot(1,3,3)
        plot(y,z)
        xlabel('y')
        ylabel('z')
        
        figure(3*n+9)
        plot3(x,y,z)
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
end

%% 24.74<r<28
for n=1:3
    r = random('Uniform',28,50,[1 1]);
    for j=1:2
        x(1) = random('Uniform',1,10,[1 1]);
        y(1) = random('Uniform',1,10,[1 1]);
        z(1) = random('Uniform',1,10,[1 1]);
        for i=1:length(t)-1
            dx = s*(y(i)-x(i))*dt;
            dy = (r*x(i)-y(i)-x(i)*z(i))*dt;
            dz = (x(i)*y(i)-b*z(i))*dt;
            
            x(i+1) = x(i) + dx;
            y(i+1) = y(i) + dy;
            z(i+1) = z(i) + dz;
        end
        figure(3*n+16)
        subplot(3,2,j)
        plot(t,x)
        title(['r=',num2str(r),',x0=',num2str(x(1)),',y0=',num2str(y(1)),',z0=',num2str(z(1))])
        xlabel('t')
        ylabel('x')
        subplot(3,2,2+j)
        plot(t,y)
        xlabel('t')
        ylabel('y')
        subplot(3,2,4+j)
        plot(t,z)
        xlabel('t')
        ylabel('z')
        
        figure(3*n+17)
        subplot(1,3,1)
        plot(x,y)
        xlabel('x')
        ylabel('y')
        subplot(1,3,2)
        plot(x,z)
        xlabel('x')
        ylabel('z')
        subplot(1,3,3)
        plot(y,z)
        xlabel('y')
        ylabel('z')
        
        figure(3*n+18)
        plot3(x,y,z)
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
end

%%

clc; clear all; close all;
%% classification

load('1_sig10_250.mat')
lyapExp(10) = lyapunovExponent(rec_signal,250);
approxEnt(10) = approximateEntropy(rec_signal);
corDim(10) = correlationDimension(rec_signal);
Label(10) = 1;

load('1_sig9_250.mat')
lyapExp(9) = lyapunovExponent(rec_signal,250);
approxEnt(9) = approximateEntropy(rec_signal);
corDim(9) = correlationDimension(rec_signal);
Label(9) = 1;

load('1_sig8_250.mat')
lyapExp(8) = lyapunovExponent(rec_signal,250);
approxEnt(8) = approximateEntropy(rec_signal);
corDim(8) = correlationDimension(rec_signal);
Label(8) = 1;

load('1_sig7_256.mat')
lyapExp(7) = lyapunovExponent(rec_signal,256);
approxEnt(7) = approximateEntropy(rec_signal);
corDim(7) = correlationDimension(rec_signal);
Label(7) = 1;

load('1_sig6_256.mat')
lyapExp(6) = lyapunovExponent(rec_signal,256);
approxEnt(6) = approximateEntropy(rec_signal);
corDim(6) = correlationDimension(rec_signal);
Label(6) = 1;

load('1_sig5_256.mat')
lyapExp(5) = lyapunovExponent(rec_signal,256);
approxEnt(5) = approximateEntropy(rec_signal);
corDim(5) = correlationDimension(rec_signal);
Label(5) = 1;

load('1_sig4_256.mat')
lyapExp(4) = lyapunovExponent(rec_signal,256);
approxEnt(4) = approximateEntropy(rec_signal);
corDim(4) = correlationDimension(rec_signal);
Label(4) = 1;

load('1_sig3_250.mat')
lyapExp(3) = lyapunovExponent(rec_signal,250);
approxEnt(3) = approximateEntropy(rec_signal);
corDim(3) = correlationDimension(rec_signal);
Label(3) = 1;

load('1_sig2_256.mat')
lyapExp(2) = lyapunovExponent(rec_signal,256);
approxEnt(2) = approximateEntropy(rec_signal);
corDim(2) = correlationDimension(rec_signal);
Label(2) = 1;

load('1_sig1_250.mat')
lyapExp(1) = lyapunovExponent(rec_signal,250);
approxEnt(1) = approximateEntropy(rec_signal);
corDim(1) = correlationDimension(rec_signal);
Label(1) = 1;

load('0_sig10_250.mat')
lyapExp(20) = lyapunovExponent(rec_signal,250);
approxEnt(20) = approximateEntropy(rec_signal);
corDim(20) = correlationDimension(rec_signal);
Label(20) = 0;

load('0_sig9_256.mat')
lyapExp(19) = lyapunovExponent(rec_signal,256);
approxEnt(19) = approximateEntropy(rec_signal);
corDim(19) = correlationDimension(rec_signal);
Label(19) = 0;

load('0_sig8_256.mat')
lyapExp(18) = lyapunovExponent(rec_signal,256);
approxEnt(18) = approximateEntropy(rec_signal);
corDim(18) = correlationDimension(rec_signal);
Label(18) = 0;

load('0_sig7_250.mat')
lyapExp(17) = lyapunovExponent(rec_signal,250);
approxEnt(17) = approximateEntropy(rec_signal);
corDim(17) = correlationDimension(rec_signal);
Label(17) = 0;

load('0_sig6_250.mat')
lyapExp(16) = lyapunovExponent(rec_signal,250);
approxEnt(16) = approximateEntropy(rec_signal);
corDim(16) = correlationDimension(rec_signal);
Label(16) = 0;

load('0_sig5_256.mat')
lyapExp(15) = lyapunovExponent(rec_signal,256);
approxEnt(15) = approximateEntropy(rec_signal);
corDim(15) = correlationDimension(rec_signal);
Label(15) = 0;

load('0_sig4_256.mat')
lyapExp(14) = lyapunovExponent(rec_signal,256);
approxEnt(14) = approximateEntropy(rec_signal);
corDim(14) = correlationDimension(rec_signal);
Label(14) = 0;

load('0_sig3_250.mat')
lyapExp(13) = lyapunovExponent(rec_signal,250);
approxEnt(13) = approximateEntropy(rec_signal);
corDim(13) = correlationDimension(rec_signal);
Label(13) = 0;

load('0_sig2_256.mat')
lyapExp(12) = lyapunovExponent(rec_signal,256);
approxEnt(12) = approximateEntropy(rec_signal);
corDim(12) = correlationDimension(rec_signal);
Label(12) = 0;

load('0_sig1_250.mat')
lyapExp(11) = lyapunovExponent(rec_signal,250);
approxEnt(11) = approximateEntropy(rec_signal);
corDim(11) = correlationDimension(rec_signal);
Label(11) = 0;

%%
class0 = find(Label==0);
class1 = find(Label==1);
Mean_class0 = [mean(lyapExp(class0)),mean(approxEnt(class0)),mean(corDim(class0))]
Mean_class1 = [mean(lyapExp(class1)),mean(approxEnt(class1)),mean(corDim(class1))]
STD_class0 = [std(lyapExp(class0)),std(approxEnt(class0)),std(corDim(class0))]
STD_class1 = [std(lyapExp(class1)),std(approxEnt(class1)),std(corDim(class1))]

N_lyapExp = (lyapExp-mean(lyapExp))/std(lyapExp);
N_approxEnt = (approxEnt-mean(approxEnt))/std(approxEnt);
N_corDim = (corDim-mean(corDim))/std(corDim);
N_Mean_class0 = [mean(N_lyapExp(class0)),mean(N_approxEnt(class0)),mean(N_corDim(class0))]
N_Mean_class1 = [mean(N_lyapExp(class1)),mean(N_approxEnt(class1)),mean(N_corDim(class1))]
N_STD_class0 = [std(N_lyapExp(class0)),std(N_approxEnt(class0)),std(N_corDim(class0))]
N_STD_class1 = [std(N_lyapExp(class1)),std(N_approxEnt(class1)),std(N_corDim(class1))]

figure(1)
scatter3(lyapExp,approxEnt,corDim,20,Label)

%% box-counting
clc; clear all; close all;

%% koch_snowflake
c = imread('koch_snowflake.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_koch_snowflake] = boxcount(c,'slope');
Fractal_koch_snowflake = s_koch_snowflake(1);
%% Sierpinski_carpet
c = imread('Sierpinski_carpet.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_Sierpinski_carpet] = boxcount(c,'slope');
Fractal_Sierpinski_carpet = s_Sierpinski_carpet(1);
%% Sierpinski
c = imread('Sierpinski.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_Sierpinski] = boxcount(c,'slope');
Fractal_Sierpinski_triangle = s_Sierpinski(1);

%% logistic-map
c = imread('logistic-map.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_logistic_map] = boxcount(c,'slope');
Fractal_logistic_map = s_logistic_map(1);

%% Cantor-set
c = imread('Cantor-set.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_Cantor_set] = boxcount(c,'slope');
Fractal_Cantor_set = s_Cantor_set(1);

%%
c = imread('1.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_1] = boxcount(c,'slope');
Fractal1 = s_1(1);

%%
c = imread('2.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_2] = boxcount(c,'slope');
Fractal2 = s_2(1);

%%
c = imread('3.jpg');
c = (c<198);
figure
imagesc(~c)
colormap gray
axis square
figure
boxcount(c)
figure
[n,r,s_3] = boxcount(c,'slope');
Fractal3 = s_3(1);

