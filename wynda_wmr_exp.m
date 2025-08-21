%% Research code by Agus Hasan & Shinta Kosmaga
clear all;
close all;

%% number of variables and coefficients
n = 5;
r = 20;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 0;

%% true parameters
m  = 0.85613;
lf = 0.06874;
lr = 0.06726;
Kf = 0.04;
Kr = 0.0435;
I = 0.00794;
l = lf+lr;

%% state initialization
x        = zeros(n,1);
xbar     = x;
y        = x;
thetabar = zeros(r,1);
Kfbar    = Kf;
Krbar    = Kr;
Ibar     = zeros(1,1);
lfbar    = zeros(1,1);
lrbar    = zeros(1,1);
ebar     = xbar-x;
t        = 0;
mbar     = 0;

%% initial control inputs
V= 0;
S= 0;
B= 0;

%% for plotting
SArray         = [];
VArray         = [];
BArray         = [];
xArray         = [];
xbarArray      = [];
yArray         = [];
thetabarArray  = [];
KfbarArray     = [];
KrbarArray     = [];
IbarArray      = [];
lfbarArray     = [];
lrbarArray     = [];
ebarArray      = [];
tArray         = [];
mbarArray      = [];
BArray         = [];

%% initialization for estimator
lambdav = 0.7;
lambdat = 0.995;
Rx      = 1e5*eye(n);
Rt      = 1e5*eye(n);
Px      = 0.01*eye(n);
Pt      = 1.35e+04*eye(r);
Gamma   = 1*zeros(n,r);

%% pre-allocate data logging
data_log = readmatrix('Data_Exp.xlsx');

%% data processing
%filter data
data_log(:,4)=sgolayfilt(data_log(:,4), 4, 151); %x
data_log(:,5)=sgolayfilt(data_log(:,5), 4, 151); %y
data_log(:,6)=sgolayfilt(data_log(:,6), 4, 151); %theta

%shifting x and y to origin (0,0)
x_raw=data_log(:,4);
y_raw=data_log(:,5);
data_log(:,4)=x_raw-x_raw(1);
data_log(:,5)=y_raw-y_raw(1);

%xdot, ydot, thdot
dx=gradient(data_log(:,4))./gradient(data_log(:,1));
dy=gradient(data_log(:,5))./gradient(data_log(:,1));
dth=gradient(data_log(:,6))./gradient(data_log(:,1));

%beta
gam=unwrap(atan2(dy,dx));
gam(gam < 0) = gam(gam < 0) + 2*pi();
beta=gam-data_log(:,6);
beta(beta < 0) = beta(beta < 0) + 2*pi();

%filter data
data_log(:,7)=sgolayfilt(dy, 4, 151); %dy
data_log(:,8)=sgolayfilt(dth, 4, 151); %dth
data_log(:,10)=sgolayfilt(beta, 4, 151); %beta

%% Loop Real-time
for i=1:1:size(data_log,1)
        
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];      
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar];
    KfbarArray     = [KfbarArray Kfbar];
    KrbarArray     = [KrbarArray Krbar];
    IbarArray      = [IbarArray Ibar];
    lfbarArray     = [lfbarArray lfbar];
    lrbarArray     = [lrbarArray lrbar];
    ebarArray      = [ebarArray ebar];
    VArray         = [VArray V];
    SArray         = [SArray S];
    tArray         = [tArray t];
    mbarArray      = [mbarArray mbar];
    BArray         = [BArray B];
    
    % real-time data
    x(1) = data_log(i,5);  % y
    x(2) = data_log(i,4);  % x
    x(3) = data_log(i,6);  % theta
    x(4) = data_log(i,7);  % y dot
    x(5) = data_log(i,8);  % theta_dot
    B = data_log(i,10);    % beta
    V = data_log(i,2);     % velocity
    S = -data_log(i,3);     % steer
    delta_t = data_log(i,9);
    t = tArray(end)+delta_t;
           
    y = C*x+delta_t*R^2*randn(n,1);

    Phi = [V*sin(B+y(3)) 0 0 0 zeros(16,1)';
        zeros(4,1)' V*cos(B+y(3)) 0 0 0 zeros(12,1)'
        zeros(8,1)' y(5) 0 0 0 zeros(8,1)'
        zeros(12,1)' y(3) y(4)/V y(5)/V S zeros(4,1)'
        zeros(16,1)' y(3) y(4)/V y(5)/V S];
         
    % estimation using adaptive observer
    Kx = Px*C'*inv(C*Px*C'+Rx);
    Kt = Pt*Gamma'*C'*inv(C*Gamma*Pt*Gamma'*C'+Rt);
    Gamma = (eye(n)-Kx*C)*Gamma;
    
    xbar = xbar+(Kx+Gamma*Kt)*(y-C*xbar);
    thetabar = thetabar-Kt*(y-C*xbar);
    
    xbar = A*xbar+Phi*thetabar;
    thetabar = thetabar;
    
    Px = (1/lambdav)*eye(n)*(eye(n)-Kx*C)*Px*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*C*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

    % estimation result
    Kfbar = Kf;
    Krbar = Kr;
    mbar  = (2*(Kfbar)*delta_t)/thetabar(16);
    lfbar = (2*Krbar*l-(thetabar(15)*m/delta_t))/(2*Kfbar+2*Krbar);
    lrbar = l-lfbar;
    Ibar  = (2*Kfbar*lfbar*delta_t)/thetabar(20);

    ebar = xbar-x;
end

%% final value and RMSE
mbarA=mean(mbarArray(end-150:end));
lfbarA=mean(lfbarArray(end-150:end));
lrbarA=mean(lrbarArray(end-150:end));
IbarA=mean(IbarArray(end-150:end));

mbarerr = sqrt(mean((mbarArray(end-150:end) - m).^2));
lfbarerr = sqrt(mean((lfbarArray(end-150:end) - lf).^2));
lrbarerr = sqrt(mean((lrbarArray(end-150:end) - lr).^2));
Ibarerr = sqrt(mean((IbarArray(end-150:end) - I).^2));

%% plotting
figure(1)
subplot(2,2,1)
plot(tArray,(m)*ones(1,length(tArray)),'-','LineWidth',10);
hold on;
plot(tArray,mbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
legend('3D Model','WyNDA')
grid on;
grid minor;
ylabel('m (kg)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-5 5])
subplot(2,2,2)
plot(tArray,I*ones(1,length(tArray)),'-','LineWidth',10);
hold on;
plot(tArray,IbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('Iy (kgm^2)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-1 1])
subplot(2,2,3)
plot(tArray,lf*ones(1,length(tArray)),'-','LineWidth',10);
hold on;
plot(tArray,lfbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('df (m)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-1 1])
subplot(2,2,4)
plot(tArray,lr*ones(1,length(tArray)),'-','LineWidth',10);
hold on;
plot(tArray,lrbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('dr (m)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-1 1])

figure(2)
subplot(2,1,1)
plot(tArray,SArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\gamma (rad)','FontSize',18)
xlabel('t (s)','FontSize',18)
subplot(2,1,2)
plot(tArray,VArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('v (m/s)','FontSize',18)
xlabel('t (s)','FontSize',18)

figure(3)
plot(xArray(2,:),xArray(1,:),'-k','LineWidth',10);
hold on;
plot(xbarArray(2,:),xbarArray(1,:),':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
legend('Eksperimen','WyNDA')
ylabel('y (m)','FontSize',18)
xlabel('x (m)','FontSize',18)
ylim([-0.8 2])
xlim([-0.8 2])