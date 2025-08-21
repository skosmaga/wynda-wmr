%% Research code by Agus Hasan & Shinta Kosmaga
clear all;
close all;

%% number of variables and coefficients
n = 5;
r = 20;

%% time horizon
tf  = 50;
dt  = 0.05;
t   = dt:dt:tf;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 0;

%% state initialization
x        = zeros(n,1);
xbar     = x;
y        = x;
thetabar = zeros(r,1);
Kfbar    = zeros(1,1);
Krbar    = zeros(1,1);
Ibar     = zeros(1,1);
lfbar    = zeros(1,1);
lrbar    = zeros(1,1);
ebar     = xbar-x;
mbar     = 0;

%% true parameters
m  = 0.85613;
lf = 0.06874;
lr = 0.06726;
Kf = 0.04;
Kr = 0.0435;
I = 0.00794;
l = lr+lf;

%% initial control inputs
S = 0;
V = 0;
B = 0;

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
mbarArray      = [];

%% initializaion for estimator
lambdav = 0.9999;
lambdat = 0.9995;
Rx      = 10*eye(n);
Rt      = 1*eye(n);
Px      = 1*eye(n);
Pt      = 1*eye(r);
Gamma   = 1*zeros(n,r);

%% simulation
for i=dt:dt:tf
        
    VArray         = [VArray V];
    SArray         = [SArray S];
    BArray         = [BArray B];
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
    mbarArray      = [mbarArray mbar];

    if i<15
        V = 0.1;
        S = -0.2;
    else
        V = 0.14;
        S = 0.2;
    end
    
    % simulate the system using Runge-Kutta
    a1=2*(Kf+Kr);
    a2=2*((lf*Kf)-(lr*Kr));
    a3=2*((lf^2)*Kf)+((lr^2)*Kr);
    a4=2*lf*Kf;

    k1=V*sin(B+x(3));
    l1=V*cos(B+x(3));
    m1=x(5);
    n1=(a1*x(3)/m)-(a1*x(4)/(m*V))-(a2*x(5)/(m*V))+(2*Kf*S/m);
    o1=(a2*x(3)/I)-(a2*x(4)/(I*V))-(a3*x(5)/(I*V))+(a4*S/I);

    k2=V*sin(B+(x(3)+(0.5*dt*m1)));
    l2=V*cos(B+(x(3)+(0.5*dt*m1)));
    m2=x(5)+(0.5*dt*o1);
    n2=(a1*(x(3)+(0.5*dt*m1))/m)-(a1*(x(4)+(0.5*dt*n1))/(m*V))-(a2*(x(5)+(0.5*dt*o1))/(m*V))+(2*Kf*S/m);
    o2=(a2*(x(3)+(0.5*dt*m1))/I)-(a2*(x(4)+(0.5*dt*n1))/(I*V))-(a3*(x(5)+(0.5*dt*o1))/(I*V))+(a4*S/I);

    k3=V*sin(B+(x(3)+(0.5*dt*m2)));
    l3=V*cos(B+(x(3)+(0.5*dt*m2)));
    m3=x(5)+(0.5*dt*o2);
    n3=(a1*(x(3)+(0.5*dt*m2))/m)-(a1*(x(4)+(0.5*dt*n2))/(m*V))-(a2*(x(5)+(0.5*dt*o2))/(m*V))+(2*Kf*S/m);
    o3=(a2*(x(3)+(0.5*dt*m2))/I)-(a2*(x(4)+(0.5*dt*n2))/(I*V))-(a3*(x(5)+(0.5*dt*o2))/(I*V))+(a4*S/I);
    
    k4=V*sin(B+(x(3)+(dt*m3)));
    l4=V*cos(B+(x(3)+(dt*m3)));
    m4=x(5)+(dt*o3);
    n4=(a1*(x(3)+(dt*m3))/m)-(a1*(x(4)+(dt*n3))/(m*V))-(a2*(x(5)+(dt*o3))/(m*V))+(2*Kf*S/m);
    o4=(a2*(x(3)+(dt*m3))/I)-(a2*(x(4)+(dt*n3))/(I*V))-(a3*(x(5)+(dt*o3))/(I*V))+(a4*S/I);
    
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    x(3) = x(3) + (dt/6)*(m1+2*m2+2*m3+m4);
    x(4) = x(4) + (dt/6)*(n1+2*n2+2*n3+n4);
    x(5) = x(5) + (dt/6)*(o1+2*o2+2*o3+o4);
    
    y = C*x+dt*R^2*randn(n,1);

    Phi = [V*sin(B+y(3)) 0 0 0 zeros(16,1)';
        zeros(4,1)' V*cos(B+y(3)) 0 0 0 zeros(12,1)'
        zeros(8,1)' y(5) 0 0 0 zeros(8,1)'
        zeros(12,1)' y(3) y(4)/V y(5)/V S zeros(4,1)'
        zeros(16,1)' y(3) y(4)/V y(5)/V S];

    % Estimation using adaptive observer
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
    mbar  = 2*Kfbar*dt/thetabar(16);
    lfbar = (2*Krbar*l-(thetabar(15)*m/dt))/(2*Kfbar+2*Krbar);
    lrbar = l-lfbar;
    Ibar  = (2*Kfbar*lfbar*dt)/thetabar(20);
        
    ebar = xbar-x;
    B=unwrap(atan2(x(4),V));
end

%% final value and RMSE
mbar_avg = mean(mbarArray(850:end));
Ibar_avg = mean(IbarArray(850:end));
lfbar_avg = mean(lfbarArray(850:end));
lrbar_avg = mean(lrbarArray(850:end));

mbarerr = sqrt(mean((mbarArray(850:end) - m).^2));
lfbarerr = sqrt(mean((lfbarArray(850:end) - lf).^2));
lrbarerr = sqrt(mean((lrbarArray(850:end) - lr).^2));
Iarerr = sqrt(mean((IbarArray(850:end) - I).^2));

%% plotting
figure(1)
subplot(2,2,1)
plot(t,(m)*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,mbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
legend('3D Model','WyNDA')
grid on;
grid minor;
ylabel('m (kg)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-3 3])
subplot(2,2,2)
plot(t,I*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,IbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('Iy (kgm^2)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-0.1 0.1])
subplot(2,2,3)
plot(t,lf*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,lfbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('df (m)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-1 1])
subplot(2,2,4)
plot(t,lr*ones(1,length(t)),'-','LineWidth',10);
hold on;
plot(t,lrbarArray,':','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('lr (m)','FontSize',18)
xlabel('t (s)','FontSize',18)
ylim([-1 1])

figure(2)
subplot(2,1,1)
plot(t,SArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
ylabel('\gamma (rad)','FontSize',18)
xlabel('t (s)','FontSize',18)
subplot(2,1,2)
plot(t,VArray,'-k','LineWidth',10);
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
ylim([-2 0.5])
xlim([-2 0.5])