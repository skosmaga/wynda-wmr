%% ============================================================
%  Research Code
%  Authors : Agus Hasan & Shinta Kosmaga
%  Purpose : WMR parameter estimation using WyNDA and SINDy
%% ============================================================
clear all; close all;

%% Number of variables and coefficients
n = 5;
r = 20;

%% Time horizon
tf  = 50;
dt  = 0.05;
t   = dt:dt:tf;
nsteps = length(t);

%% System description
A = eye(n);
C = eye(n);

%% Noise
R = 0.15;

%% State initialization
x        = zeros(n,1);
xbar     = x;
y        = x;
thetabar = zeros(r,1);
Ibar     = zeros(1,1);
dfbar    = zeros(1,1);
drbar    = zeros(1,1);
mbar     = zeros(1,1);
ebar     = xbar-x;
mbarS     = 0;
IbarS     = zeros(1,1);
dfbarS    = zeros(1,1);
drbarS    = zeros(1,1);

%% True WMR parameters
m  = 0.85613;
df = 0.06874;
dr = 0.06726;
Kf = 0.04;
Kr = 0.0435;
I = 0.00794;
d = dr+df;

%% Inputs initialization
S = 0;
V = 0;
B = 0;

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
mbarSArray     = [];
IbarSArray     = [];
dfbarSArray    = [];
drbarSArray    = [];
dy             = [];
Xi_hist        = [];
t_hist         = [];

%% WyNDA Hyperparameter Initialization
lambdav = 0.9999;
lambdat = 0.9995;
Rx      = 10*eye(n);
Rt      = 1*eye(n);
Px      = 1*eye(n);
Pt      = 1*eye(r);
Gamma   = 1*zeros(n,r);

%% WyNDA simulation loop
for i=1:nsteps

    % Stored data
    VArray         = [VArray V];
    SArray         = [SArray S];
    BArray         = [BArray B];
    xArray         = [xArray x];   
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    xbarArray      = [xbarArray xbar];  
    IbarArray      = [IbarArray Ibar];
    dfbarArray     = [dfbarArray dfbar];
    drbarArray     = [drbarArray drbar];
    mbarArray      = [mbarArray mbar];
    ebarArray      = [ebarArray ebar];

    % Input profile
    if i<300
        V = 0.1;
        S = -0.2;
    else
        V = 0.14;
        S = 0.2;
    end

    % Simulate the system using Runge-Kutta
    a1=2*(Kf+Kr);
    a2=2*((df*Kf)-(dr*Kr));
    a3=2*((df^2)*Kf)+((dr^2)*Kr);
    a4=2*df*Kf;

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

    % Measurement data    
    y = C*x+dt*R^2*randn(n,1);
    
    % Approximation function
    Phi = [V*sin(B+y(3)) 0 0 0 zeros(16,1)';
        zeros(4,1)' V*cos(B+y(3)) 0 0 0 zeros(12,1)'
        zeros(8,1)' y(5) 0 0 0 zeros(8,1)'
        zeros(12,1)' y(3) y(4)/V y(5)/V S zeros(4,1)'
        zeros(16,1)' y(3) y(4)/V y(5)/V S];
    
    % True model dynamics (for SINDy)
    f(:,i) = [V*sin(B+y(3));
              V*cos(B+y(3));
              y(5);
              0;
              0];

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

    % Parameter Recovery
    mbar  = 2*Kf*dt/thetabar(16);
    dfbar = (2*Kr*d-(thetabar(15)*mbar/dt))/(2*Kf+2*Kr);
    drbar = d-dfbar;
    Ibar  = (2*Kf*dfbar*dt)/thetabar(20);

    ebar = xbar-x;
    B=unwrap(atan2(x(4),V));
end

%% SINDy simulation
% Numerical differentiation
for i = 3:length(yArray)-3
    dy(:,i-2) = (1/(12*dt)) * (-yArray(:,i+2) + 8*yArray(:,i+1) ...
                               - 8*yArray(:,i-1) + yArray(:,i-2));
end
res = dy - f(:,3:end-3);  % Residuals for parameter estimation
yo = yArray(:,3:end-3);   % Aligned measurements
Vo = VArray(:,3:end-3);   % Aligned inputs
So = SArray(:,3:end-3);

% Library construction
Theta = [ yo(3,:);
          yo(4,:)./Vo;
          yo(5,:)./Vo;
          So]';

% Mask definition
% Target states: x4, x5
mask = false(4,2);
mask(1:4,1:2) = true;

lambda = 0;  % No thresholding
nDelay=size(yArray,2)-size(Theta,1);

% Initialization
for k = 1:1:nDelay
    KfbarSArray(k)=0;
    KrbarSArray(k)=0;
    mbarSArray(k)=0;
    dfbarSArray(k)=0;
    drbarSArray(k)=0;
    IbarSArray(k)=0;
end


% Sparse Regression
for k = nDelay+1:1:size(yArray,2)
    IbarSArray      = [IbarSArray IbarS];
    dfbarSArray     = [dfbarSArray dfbarS];
    drbarSArray     = [drbarSArray drbarS];
    mbarSArray      = [mbarSArray mbarS];

    Theta_k = Theta(1:k-nDelay,:);
    res_k   = res(:,1:k-nDelay)';   % (k x nStates)

    Xi_k = sparsifyDynamics(Theta_k, res_k(:,4:5), lambda, mask);

    Xi_hist(:,end+1) = [Xi_k(1:4,1)',Xi_k(1:4,2)'];
    t_hist(end+1) = k*dt;

    mbarS  = 2*Kf/Xi_hist(4,end);
    dfbarS = (2*Kr*d-(Xi_hist(3,end)*mbarS))/(2*Kf+2*Kr);
    drbarS = d-dfbarS;
    IbarS  = (2*Kf*dfbarS)/Xi_hist(8,end);
end


%% Final Estimation Statistics (Steady-State)

mbar_avg = mean(mbarArray(850:end));
Ibar_avg = mean(IbarArray(850:end));
dfbar_avg = mean(dfbarArray(850:end));
drbar_avg = mean(drbarArray(850:end));

mbarerr = sqrt(mean((mbarArray(850:end) - m).^2));
dfbarerr = sqrt(mean((dfbarArray(850:end) - df).^2));
drbarerr = sqrt(mean((drbarArray(850:end) - dr).^2));
Iarerr = sqrt(mean((IbarArray(850:end) - I).^2));

mbarerrp = mean(abs((m-mbarArray(850:end))./m)) * 100;
dfbarerrp = mean(abs((df-dfbarArray(850:end))./df)) * 100;
drbarerrp = mean(abs((dr-drbarArray(850:end))./dr)) * 100;
Ibarerrp = mean(abs((I-IbarArray(850:end))./I)) * 100;

mbarS_avg = mean(mbarSArray(850:end));
IbarS_avg = mean(IbarSArray(850:end));
dfbarS_avg = mean(dfbarSArray(850:end));
drbarS_avg = mean(drbarSArray(850:end));

mbarSerr = sqrt(mean((mbarSArray(850:end) - m).^2));
dfbarSerr = sqrt(mean((dfbarSArray(850:end) - df).^2));
drbarSerr = sqrt(mean((drbarSArray(850:end) - dr).^2));
IarSerr = sqrt(mean((IbarSArray(850:end) - I).^2));

mbarSerrp = mean(abs((m-mbarSArray(850:end))./m)) * 100;
dfbarSerrp = mean(abs((df-dfbarSArray(850:end))./df)) * 100;
drbarSerrp = mean(abs((dr-drbarSArray(850:end))./dr)) * 100;
IbarSerrp = mean(abs((I-IbarSArray(850:end))./I)) * 100;

%% Plotting

% Parameter estimation (WyNDA)
figure(1)

paramTrue = {m, I, df, dr};
paramWy   = {mbarArray, IbarArray, dfbarArray, drbarArray};
ylab      = {'m (kg)','I_y (kg m^2)','d_f (m)','d_r (m)'};
ylims     = {[-3 3], [-0.1 0.1], [-1 1], [-1 1]};

for i = 1:4
    subplot(2,2,i)
    plot(t,paramTrue{i}*ones(1,nsteps),'-','LineWidth',10); hold on;
    plot(t,paramWy{i},':','LineWidth',10);
    ylim(ylims{i})
    formatAxis('t (s)', ylab{i})
end

% Parameter estimation (SINDy)
figure(2)

paramTrue = {m, I, df, dr};
paramSi   = {mbarSArray, IbarSArray, dfbarSArray, drbarSArray};
ylab      = {'m (kg)','I_y (kg m^2)','d_f (m)','d_r (m)'};
ylims     = {[-3 3], [-0.1 0.1], [-1 1], [-1 1]};

for i = 1:4
    subplot(2,2,i)
    plot(t,paramTrue{i}*ones(1,nsteps),'-','LineWidth',10); hold on;
    plot(t,paramSi{i},':','LineWidth',10);
    ylim(ylims{i})
    formatAxis('t (s)', ylab{i})
end

% Input
figure(3)
subplot(2,1,1)
plot(t,SArray,'-k','LineWidth',10);
formatAxis('t (s)','\gamma (rad)')

subplot(2,1,2)
plot(t,VArray,'-k','LineWidth',10);
formatAxis('t (s)','v (m/s)')

figure(5)
plot(xArray(2,:),xArray(1,:),'-k','LineWidth',10); hold on;
plot(xbarArray(2,:),xbarArray(1,:),':','LineWidth',10);
legend('Eksperimen','WyNDA')
formatAxis('x (m)','y (m)')
xlim([-2 0.5]); ylim([-2 0.5])

function formatAxis(xlab, ylab)
    set(gca,'Color','white','LineWidth',3,'FontSize',12)
    grid on; grid minor;
    xlabel(xlab,'FontSize',18)
    ylabel(ylab,'FontSize',18)
end
