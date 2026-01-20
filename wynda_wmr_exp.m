%% ============================================================
%  Research Code
%  Authors : Agus Hasan & Shinta Kosmaga
%  Purpose : WMR parameter estimation using WyNDA (Experimental)
%% ============================================================
clear all; close all;

%% Number of variables and coefficients
n = 5;
r = 20;

%% System description
A = eye(n);
C = eye(n);

%% Noise
R = 0;

%% True WMR parameters
m  = 0.85613;
df = 0.06874;
dr = 0.06726;
Kf = 0.04;
Kr = 0.0435;
I = 0.00794;
d = df+dr;

%% State initialization
x        = zeros(n,1);
xbar     = x;
y        = x;
thetabar = zeros(r,1);
Ibar     = zeros(1,1);
dfbar    = zeros(1,1);
drbar    = zeros(1,1);
mbar     = 0;
ebar     = xbar-x;
t        = 0;

%% Inputs initialization
V= 0;
S= 0;
B= 0;

%% Data storage
SArray         = [];
VArray         = [];
BArray         = [];
xArray         = [];
xbarArray      = [];
yArray         = [];
thetabarArray  = [];
IbarArray      = [];
dfbarArray     = [];
drbarArray     = [];
mbarArray      = [];
ebarArray      = [];
tArray         = [];
BArray         = [];

%% WyNDA Hyperparameter Initialization
lambdav = 0.7;
lambdat = 0.995;
Rx      = 1e5*eye(n);
Rt      = 1e5*eye(n);
Px      = 0.01*eye(n);
Pt      = 1.35e+04*eye(r);
Gamma   = 1*zeros(n,r);

%% Pre-allocate for data logging
data_log = readmatrix('Cek8.xlsx');

%% Data processing
% Filter data
data_log(:,4)=sgolayfilt(data_log(:,4), 4, 151); %x
data_log(:,5)=sgolayfilt(data_log(:,5), 4, 151); %y
data_log(:,6)=sgolayfilt(data_log(:,6), 4, 151); %theta

% Shifting x and y to origin (0,0)
x_raw=data_log(:,4);
y_raw=data_log(:,5);
data_log(:,4)=x_raw-x_raw(1);
data_log(:,5)=y_raw-y_raw(1);

% xdot, ydot, thdot
dx=gradient(data_log(:,4))./gradient(data_log(:,1));
dy=gradient(data_log(:,5))./gradient(data_log(:,1));
dth=gradient(data_log(:,6))./gradient(data_log(:,1));

% beta
gam=unwrap(atan2(dy,dx));
gam(gam < 0) = gam(gam < 0) + 2*pi();
beta=gam-data_log(:,6);
beta(beta < 0) = beta(beta < 0) + 2*pi();

% Filter data
data_log(:,7)=sgolayfilt(dy, 4, 151); %dy
data_log(:,8)=sgolayfilt(dth, 4, 151); %dth
data_log(:,10)=sgolayfilt(beta, 4, 151); %beta

%% WyNDA simulation loop
for i=1:1:size(data_log,1)
    
    % Stored data
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];      
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar];
    IbarArray      = [IbarArray Ibar];
    dfbarArray     = [dfbarArray dfbar];
    drbarArray     = [drbarArray drbar];
    mbarArray      = [mbarArray mbar];
    ebarArray      = [ebarArray ebar];
    VArray         = [VArray V];
    SArray         = [SArray S];
    tArray         = [tArray t];
    BArray         = [BArray B];
    
    % Real-time data
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
    
    % Measurement data 
    y = C*x+delta_t*R^2*randn(n,1);

    % Approximation function
    Phi = [V*sin(B+y(3)) 0 0 0 zeros(16,1)';
        zeros(4,1)' V*cos(B+y(3)) 0 0 0 zeros(12,1)'
        zeros(8,1)' y(5) 0 0 0 zeros(8,1)'
        zeros(12,1)' y(3) y(4)/V y(5)/V S zeros(4,1)'
        zeros(16,1)' y(3) y(4)/V y(5)/V S];
         
    % Adaptive observer
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
    mbar  = (2*(Kf)*delta_t)/thetabar(16);
    dfbar = (2*Kr*d-(thetabar(15)*mbar/delta_t))/(2*Kf+2*Kr);
    drbar = d-dfbar;
    Ibar  = (2*Kf*dfbar*delta_t)/thetabar(20);

    ebar = xbar-x;
end

%% Final Estimation Statistics (Steady-State)

mbar_avg=mean(mbarArray(end-150:end));
dfbar_avg=mean(dfbarArray(end-150:end));
drbar_avg=mean(drbarArray(end-150:end));
Ibar_avg=mean(IbarArray(end-150:end));

mbarerr = sqrt(mean((mbarArray(end-150:end) - m).^2));
dfbarerr = sqrt(mean((dfbarArray(end-150:end) - df).^2));
drbarerr = sqrt(mean((drbarArray(end-150:end) - dr).^2));
Ibarerr = sqrt(mean((IbarArray(end-150:end) - I).^2));

%% Plotting

% Parameter estimation
figure(1)

paramTrue = {m, I, df, dr};
paramWy   = {mbarArray, IbarArray, dfbarArray, drbarArray};
ylab      = {'m (kg)','I_y (kg m^2)','d_f (m)','d_r (m)'};
ylims     = {[-5 5], [-1 1], [-1 1], [-1 1]};

for i = 1:4
    subplot(2,2,i)
    plot(tArray,paramTrue{i}*ones(1,length(tArray)),'-','LineWidth',10); hold on;
    plot(tArray,paramWy{i},':','LineWidth',10);
    ylim(ylims{i})
    formatAxis('t (s)', ylab{i})
end

% Input
figure(2)
subplot(2,1,1)
plot(tArray,SArray,'-k','LineWidth',10);
formatAxis('t (s)','\gamma (rad)')

subplot(2,1,2)
plot(tArray,VArray,'-k','LineWidth',10);
formatAxis('t (s)','v (m/s)')

figure(3)
plot(xArray(2,:),xArray(1,:),'-k','LineWidth',10); hold on;
plot(xbarArray(2,:),xbarArray(1,:),':','LineWidth',10);
legend('Eksperimen','WyNDA')
formatAxis('x (m)','y (m)')
xlim([-0.8 2]); ylim([-0.8 2])

function formatAxis(xlab, ylab)
    set(gca,'Color','white','LineWidth',3,'FontSize',12)
    grid on; grid minor;
    xlabel(xlab,'FontSize',18)
    ylabel(ylab,'FontSize',18)
end
