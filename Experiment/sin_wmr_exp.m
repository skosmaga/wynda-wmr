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
ddy            = [];
Xi_hist        = [];
t_hist         = [];

%% Pre-allocate for data logging
data_log = readmatrix('data_exp.xlsx');

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

nsteps = size(data_log,1);

%% Loop Real-time
for i=1:1:nsteps
     
    % Stored data
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];      
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar];
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
    
    % True model dynamics
    f(:,i) = [V*sin(B+y(3));
              V*cos(B+y(3));
              y(5);
              0;
              0];
end

%% SINDy simulation
% Numerical differentiation
for i = 3:length(yArray)-3
    dt = sum(data_log(i-2:i,9))/3;
    ddy(:,i-2) = (1/(12*dt)) * (-yArray(:,i+2) + 8*yArray(:,i+1) ...
                               - 8*yArray(:,i-1) + yArray(:,i-2));
end
res = ddy - f(:,3:end-3);  % Residuals for parameter estimation
yo = yArray(:,3:end-3);   % Aligned measurements
Vo = VArray(:,3:end-3);   % Aligned inputs
So = SArray(:,3:end-3);
t_hist = tArray(:,3:end-3);

% Library construction
Theta = [ yo(3,:);
          yo(4,:)./Vo;
          yo(5,:)./Vo;
          So]';

% Mask definition
% Target states: x4, x5
mask = false(4,2);
mask(1:4,4:5) = true; 

% Sparse regiression
lambda = 0;  % No thresholding

for k = 1:1:size(Theta,1)
    IbarArray      = [IbarArray Ibar];
    dfbarArray     = [dfbarArray dfbar];
    drbarArray     = [drbarArray drbar];
    mbarArray      = [mbarArray mbar];

    Theta_k = Theta(1:k,:);
    res_k   = res(:,1:k)';   % (k x nStates)

    Xi_k = sparsifyDynamics(Theta_k, res_k, lambda, mask);

    Xi_hist(:,end+1) = [Xi_k(1:4,4)',Xi_k(1:4,5)'];

    mbar  = 2*Kf/Xi_hist(4,end);
    dfbar = (2*Kr*d-(Xi_hist(3,end)*mbar))/(2*Kf+2*Kr);
    drbar = d-dfbar;
    Ibar  = (2*Kf*dfbar)/Xi_hist(8,end);
end

%% Final Estimation Statistics (Steady-State)

mbar_avg=mean(mbarArray(end-150:end));
dfbar_avg=mean(dfbarArray(end-150:end));
drbar_avg=mean(drbarArray(end-150:end));
Ibar_avg=mean(IbarArray(end-150:end));

mbarerr = sqrt(mean((mbarArray(end-150:end) - m).^2));
Ibarerr = sqrt(mean((IbarArray(end-150:end) - I).^2));
dfbarerr = sqrt(mean((dfbarArray(end-150:end) - df).^2));
drbarerr = sqrt(mean((drbarArray(end-150:end) - dr).^2));

%% Plotting
% Parameter estimation
figure(1)

paramTrue = {m, I, df, dr};
paramWy   = {mbarArray, IbarArray, dfbarArray, drbarArray};
ylab      = {'m (kg)','I_y (kg m^2)','d_f (m)','d_r (m)'};
ylims     = {[-5 5], [-1 1], [-1 1], [-1 1]};

for i = 1:4
    subplot(2,2,i)
    plot(t_hist,paramTrue{i}*ones(1,length(t_hist)),'-','LineWidth',10); hold on;
    plot(t_hist,paramWy{i},':','LineWidth',10);
    ylim(ylims{i})
    formatAxis('t (s)', ylab{i})
end

function formatAxis(xlab, ylab)
    set(gca,'Color','white','LineWidth',3,'FontSize',12)
    grid on; grid minor;
    xlabel(xlab,'FontSize',18)
    ylabel(ylab,'FontSize',18)
end

% Custom SINDy with Structure Constraints
function Xi = sparsifyDynamics(Theta, dXdt, lambda, mask)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz

    % mask: [nTerms x nStates] logical matrix (true = parameter active)
    
    Xi = Theta \ dXdt;  % Initial least-squares guess
    
    for k = 1:10
        smallinds = abs(Xi) < lambda;
        Xi(smallinds) = 0;
        
        for i = 1:size(dXdt, 2)
            activeTerms = find(mask(:, i));
            if isempty(activeTerms)
                Xi(:, i) = 0;
            else
                Xi(:, i) = 0;
                Xi(activeTerms, i) = Theta(:, activeTerms) \ dXdt(:, i);
            end
        end
    end
end
