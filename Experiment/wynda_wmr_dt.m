%% ============================================================
%  Research Code
%  Authors : Agus Hasan & Shinta Kosmaga
%  Purpose : Digital twin of WMR
%% ============================================================
clear all; close all;

%% Number of variables and coefficients
n = 5;
r = 20;

%% Parameters

% WyNDA estimated
par_est.m  = 0.7642;
par_est.lf = 0.0663;
par_est.lr = 0.0697;
par_est.I = 0.1009;

% Model
par_mod.m  = 0.85613;
par_mod.lf = 0.06874;
par_mod.lr = 0.06726;
par_mod.I = 0.00794;

Kf = 0.04;
Kr = 0.0435;

par_est  = computeAuxParam(par_est, Kf, Kr);
par_mod = computeAuxParam(par_mod, Kf, Kr);

%% State initialization
x        = zeros(n,1);
x_m      = zeros(n,1);
err      = zeros(n,1);
errp     = zeros(n,1);
err_m    = zeros(n,1);
errp_m   = zeros(n,1);

%% Pre-allocate for data logging
data_log = readmatrix('data_dt.xlsx');

%% Data Processing
% Filter data
data_log(:,3)=sgolayfilt(data_log(:,3), 4, 21); %steering
data_log(:,4)=sgolayfilt(data_log(:,4), 4, 151); %x
data_log(:,5)=sgolayfilt(data_log(:,5), 4, 151); %y
data_log(:,6)=sgolayfilt(data_log(:,6), 4, 151); %theta

% Shifting x and y to origin (0,0)
x_raw=data_log(:,4);
y_raw=data_log(:,5);
data_log(:,4)=x_raw-x_raw(1);
data_log(:,5)=y_raw-y_raw(1);

dt  = data_log(:,9);
t   = [0; cumsum(dt(1:end-1))];

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

% Prellocation Array
N = size(data_log,1);

x_log   = [];
xm_log  = [];
err_log = [];
errp_log = [];
errm_log= [];
errpm_log= [];

V_log = zeros(1,N);
S_log = zeros(1,N);
B_log = zeros(1,N);

xt = data_log(:,[5 4 6 7 8])';

%% =======================
%  Simulation loop
% =======================
for k = 1:N

    V = data_log(k,2);
    S = -data_log(k,3);
    B = data_log(k,10);

    V_log(k)=V; S_log(k)=S; B_log(k)=B;

    if isMeasurementAvailable(k)
        x   = xt(:,k);
        x_m = xt(:,k);
    else
        x   = rk4Step(x, V, S, B, dt(k), par_est, Kf);
        x_m = rk4Step(x_m, V, S, B, dt(k), par_mod, Kf);

        err = x-xt(:,k);
        errp = err./xt(:,k);
        err_m = x_m-xt(:,k);
        errp_m = err_m./xt(:,k);

        err_log = [err_log err];
        errp_log = [errp_log errp];
        errm_log = [errm_log err_m];
        errpm_log = [errpm_log errp_m];
    end

    x_log(:,k)    = x;
    xm_log(:,k)   = x_m;
end
%% Statistics
rmse = sqrt(mean(err_log.^2, 2));
rmse_m = sqrt(mean(errm_log.^2, 2));

mape   = mean(abs(errp_log),2)*100;
mape_m = mean(abs(errpm_log),2)*100;

%% Plotting
figure(1)
subplot(3,1,1); plot(t,S_log,'k','LineWidth',4)
ylabel('\delta (rad)'); grid on

subplot(3,1,2); plot(t,V_log,'k','LineWidth',4)
ylabel('V (m/s)'); grid on

subplot(3,1,3); plot(t,B_log,'k','LineWidth',4)
ylabel('\beta (rad)'); xlabel('t (s)'); grid on
set(findall(gcf,'Type','axes'),'FontSize',14,'LineWidth',2)

figure(2)
plot(x_log(2,:),x_log(1,:),':','LineWidth',7); hold on
plot(data_log(:,4),data_log(:,5),'k','LineWidth',7)
plot(xm_log(2,:),xm_log(1,:),'-.','LineWidth',7)
legend('WyNDA','Experiment','3D Model')
xlabel('x (m)','FontSize',24); ylabel('y (m)','FontSize',24)
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on; grid minor;
ylim([-1 1.5])
xlim([-0.5 2])


function par = computeAuxParam(par, Kf, Kr)
    par.a1 = 2*(Kf+Kr);
    par.a2 = 2*((par.lf*Kf)-(par.lr*Kr));
    par.a3 = 2*((par.lf^2)*Kf + (par.lr^2)*Kr);
    par.a4 = 2*par.lf*Kf;
end

function x = rk4Step(x,V,S,B,dt,par,Kf)
    k1=V*sin(B+x(3));
    l1=V*cos(B+x(3));
    m1=x(5);
    n1=(par.a1*x(3)/par.m)-(par.a1*x(4)/(par.m*V))-(par.a2*x(5)/(par.m*V))+(2*Kf*S/par.m);
    o1=(par.a2*x(3)/par.I)-(par.a2*x(4)/(par.I*V))-(par.a3*x(5)/(par.I*V))+(par.a4*S/par.I);

    k2=V*sin(B+(x(3)+(0.5*dt*m1)));
    l2=V*cos(B+(x(3)+(0.5*dt*m1)));
    m2=x(5)+(0.5*dt*o1);
    n2=(par.a1*(x(3)+(0.5*dt*m1))/par.m)-(par.a1*(x(4)+(0.5*dt*n1))/(par.m*V))-(par.a2*(x(5)+(0.5*dt*o1))/(par.m*V))+(2*Kf*S/par.m);
    o2=(par.a2*(x(3)+(0.5*dt*m1))/par.I)-(par.a2*(x(4)+(0.5*dt*n1))/(par.I*V))-(par.a3*(x(5)+(0.5*dt*o1))/(par.I*V))+(par.a4*S/par.I);

    k3=V*sin(B+(x(3)+(0.5*dt*m2)));
    l3=V*cos(B+(x(3)+(0.5*dt*m2)));
    m3=x(5)+(0.5*dt*o2);
    n3=(par.a1*(x(3)+(0.5*dt*m2))/par.m)-(par.a1*(x(4)+(0.5*dt*n2))/(par.m*V))-(par.a2*(x(5)+(0.5*dt*o2))/(par.m*V))+(2*Kf*S/par.m);
    o3=(par.a2*(x(3)+(0.5*dt*m2))/par.I)-(par.a2*(x(4)+(0.5*dt*n2))/(par.I*V))-(par.a3*(x(5)+(0.5*dt*o2))/(par.I*V))+(par.a4*S/par.I);
    
    k4=V*sin(B+(x(3)+(dt*m3)));
    l4=V*cos(B+(x(3)+(dt*m3)));
    m4=x(5)+(dt*o3);
    n4=(par.a1*(x(3)+(dt*m3))/par.m)-(par.a1*(x(4)+(dt*n3))/(par.m*V))-(par.a2*(x(5)+(dt*o3))/(par.m*V))+(2*Kf*S/par.m);
    o4=(par.a2*(x(3)+(dt*m3))/par.I)-(par.a2*(x(4)+(dt*n3))/(par.I*V))-(par.a3*(x(5)+(dt*o3))/(par.I*V))+(par.a4*S/par.I);
    
    x(1) = x(1) + (dt/6)*(k1+2*k2+2*k3+k4);
    x(2) = x(2) + (dt/6)*(l1+2*l2+2*l3+l4);
    x(3) = x(3) + (dt/6)*(m1+2*m2+2*m3+m4);
    x(4) = x(4) + (dt/6)*(n1+2*n2+2*n3+n4);
    x(5) = x(5) + (dt/6)*(o1+2*o2+2*o3+o4);
end

function flag = isMeasurementAvailable(k)
    flag = (k<250 || (k>300 && k<500) || (k>550 && k<750) || ...
            (k>800 && k<1000));
end
