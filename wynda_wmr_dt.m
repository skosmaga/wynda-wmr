%% Research code by Agus Hasan
clear all;
close all;

%% number of variables and coefficients
n = 5;
r = 20;

%% parameters
% parameter from 3d model
m_m  = 0.85613;
lf_m = 0.06874;
lr_m = 0.06726;
I_m = 0.00794;
Kf = 0.04;
Kr = 0.0435;

a1_m=2*(Kf+Kr);
a2_m=2*((lf_m*Kf)-(lr_m*Kr));
a3_m=2*((lf_m^2)*Kf)+((lr_m^2)*Kr);
a4_m=2*lf_m*Kf;

% parameter from wynda
m  = 0.7642;
lf = 0.0645;
lr = 0.0715;
I = 0.099;

a1=2*(Kf+Kr);
a2=2*((lf*Kf)-(lr*Kr));
a3=2*((lf^2)*Kf)+((lr^2)*Kr);
a4=2*lf*Kf;

%% state initialization
x        = zeros(n,1);
x_m      = zeros(n,1);
x_ws     = zeros(n,1);
err      = zeros(n,1);
err_m    = zeros(n,1);
errp     = zeros(n,1);
errp_m   = zeros(n,1);

%% initial control inputs
V= 0;
S= 0;
B= 0;

%% for plotting
SArray         = [];
VArray         = [];
BArray         = [];
xArray         = [];
x_mArray       = [];
ebarArray      = [];
errArray       = [];
err_mArray     = [];
errpArray      = [];
errp_mArray    = [];
x_wsArray      = [];

%% pre-allocate data logging
data_log = readmatrix('Data_DT.xlsx');  % kosong, nanti diisi [timestamp, speed, nodemcu, x, y, z, angle]

%% data processing
% filter data
data_log(:,3)=sgolayfilt(data_log(:,3), 4, 21); %steering
data_log(:,4)=sgolayfilt(data_log(:,4), 4, 151); %x
data_log(:,5)=sgolayfilt(data_log(:,5), 4, 151); %y
data_log(:,6)=sgolayfilt(data_log(:,6), 4, 151); %theta


% shifting x and y to origin (0,0)
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

% filter data
data_log(:,7)=sgolayfilt(dy, 4, 151); %dy
data_log(:,8)=sgolayfilt(dth, 4, 151); %dth
data_log(:,10)=sgolayfilt(beta, 4, 151); %beta

% simulation time
dt_array = data_log(:,9);
t_input = [0; cumsum(dt_array(1:end-1))];

% initialization input
V_input = data_log(:,2);
S_input = data_log(:,3);
B_input = data_log(:,10);

% initialization state for simulation using parameter wynda
x(3) = data_log(1,6);
x(4) = data_log(1,7);
x(5) = data_log(1,8);

% initialization state for simulation using parameter 3d model
x_m(3) = data_log(1,6);
x_m(4) = data_log(1,7);
x_m(5) = data_log(1,8);

% experiment data
xt(:,1) = data_log(:,5);
xt(:,2) = data_log(:,4);
xt(:,3) = data_log(:,6);
xt(:,4) = data_log(:,7);
xt(:,5) = data_log(:,8);
xt=xt';

k=1;
%% simulation
for i=1:1:size(data_log,1)
    VArray         = [VArray V];
    BArray         = [BArray B];
    SArray         = [SArray S];
    xArray         = [xArray x];
    x_mArray       = [x_mArray x_m];
    errArray = [errArray err];
    x_wsArray = [x_wsArray x_ws];
    err_mArray = [err_mArray err_m];
    errpArray = [errpArray errp];
    errp_mArray = [errp_mArray errp_m];
    errp_mArray = [errp_mArray errp_m];
    
    % input for simulation using parameter wynda
    B = B_input(k);
    V = V_input(k);     
    S = -S_input(k);     
    dt=dt_array(k);
    
    % input for simulation using parameter 3d model
    B_m = B_input(k);
    V_m = V_input(k);
    S_m = -S_input(k);
    dt_m=dt_array(k);

    % state obtain from experiment data
    if k<250 || (k>300 && k<500) || (k>550 && k<750) || (k>800 && k<1000) || (k>1050 && k<1250) || (k>1300 && k<1500) || k>1550
        % state for simulation using parameter wynda
        x(1) = data_log(k,5);
        x(2) = data_log(k,4);
        x(3) = data_log(k,6);
        x(4) = data_log(k,7);
        x(5) = data_log(k,8);
        
        % state for simulation using parameter 3d model
        x_m(1) = data_log(k,5);
        x_m(2) = data_log(k,4);
        x_m(3) = data_log(k,6);
        x_m(4) = data_log(k,7);
        x_m(5) = data_log(k,8);
    else

    % data loss simulation using parameter 3d model
    k1_m=V_m*sin(B_m+x_m(3));
    l1_m=V_m*cos(B_m+x_m(3));
    m1_m=x_m(5);
    n1_m=(a1_m*x_m(3)/m_m)-(a1_m*x_m(4)/(m_m*V_m))-(a2_m*x_m(5)/(m_m*V_m))+(2*Kf*S_m/m_m);
    o1_m=(a2_m*x_m(3)/I_m)-(a2_m*x_m(4)/(I_m*V_m))-(a3_m*x_m(5)/(I_m*V_m))+(a4_m*S_m/I_m);

    k2_m=V_m*sin(B_m+(x_m(3)+(0.5*dt_m*m1_m)));
    l2_m=V_m*cos(B_m+(x_m(3)+(0.5*dt_m*m1_m)));
    m2_m=x_m(5)+(0.5*dt_m*o1_m);
    n2_m=(a1_m*(x_m(3)+(0.5*dt_m*m1_m))/m_m)-(a1_m*(x_m(4)+(0.5*dt_m*n1_m))/(m_m*V_m))-(a2_m*(x_m(5)+(0.5*dt_m*o1_m))/(m_m*V_m))+(2*Kf*S_m/m_m);
    o2_m=(a2_m*(x_m(3)+(0.5*dt_m*m1_m))/I_m)-(a2_m*(x_m(4)+(0.5*dt_m*n1_m))/(I_m*V_m))-(a3_m*(x_m(5)+(0.5*dt_m*o1_m))/(I_m*V_m))+(a4_m*S_m/I_m);

    k3_m=V_m*sin(B_m+(x_m(3)+(0.5*dt_m*m2_m)));
    l3_m=V_m*cos(B_m+(x_m(3)+(0.5*dt_m*m2_m)));
    m3_m=x_m(5)+(0.5*dt_m*o2_m);
    n3_m=(a1_m*(x_m(3)+(0.5*dt_m*m2_m))/m_m)-(a1_m*(x_m(4)+(0.5*dt_m*n2_m))/(m_m*V_m))-(a2_m*(x_m(5)+(0.5*dt_m*o2_m))/(m_m*V_m))+(2*Kf*S_m/m_m);
    o3_m=(a2_m*(x_m(3)+(0.5*dt_m*m2_m))/I_m)-(a2_m*(x_m(4)+(0.5*dt_m*n2_m))/(I_m*V_m))-(a3_m*(x_m(5)+(0.5*dt_m*o2_m))/(I_m*V_m))+(a4_m*S_m/I_m);
    
    k4_m=V_m*sin(B_m+(x_m(3)+(dt_m*m3_m)));
    l4_m=V_m*cos(B_m+(x_m(3)+(dt_m*m3_m)));
    m4_m=x_m(5)+(dt_m*o3_m);
    n4_m=(a1_m*(x_m(3)+(dt_m*m3_m))/m_m)-(a1_m*(x_m(4)+(dt_m*n3_m))/(m_m*V_m))-(a2_m*(x_m(5)+(dt_m*o3_m))/(m_m*V_m))+(2*Kf*S_m/m_m);
    o4_m=(a2_m*(x_m(3)+(dt_m*m3_m))/I_m)-(a2_m*(x_m(4)+(dt_m*n3_m))/(I_m*V_m))-(a3_m*(x_m(5)+(dt_m*o3_m))/(I_m*V_m))+(a4_m*S_m/I_m);
    
    x_m(1) = x_m(1) + (dt_m/6)*(k1_m+2*k2_m+2*k3_m+k4_m);
    x_m(2) = x_m(2) + (dt_m/6)*(l1_m+2*l2_m+2*l3_m+l4_m);
    x_m(3) = x_m(3) + (dt_m/6)*(m1_m+2*m2_m+2*m3_m+m4_m);
    x_m(4) = x_m(4) + (dt_m/6)*(n1_m+2*n2_m+2*n3_m+n4_m);
    x_m(5) = x_m(5) + (dt_m/6)*(o1_m+2*o2_m+2*o3_m+o4_m);


    % data loss simulation using parameter wynda
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
    
    % error
    err = x-xt(:,k);
    errp = err./xt(:,k);

    err_m = x_m-xt(:,k);
    errp_m = err_m./xt(:,k);
    end
    k=k+1;
end

% rmse and mape
rmse = sqrt(mean(errArray.^2, 2));
rmse_m = sqrt(mean(err_mArray.^2, 2));

mape = mean(abs(errpArray(:,5:end)), 2).*100;
mape_m = mean(abs(errp_mArray(:,5:end)), 2).*100;

%% plotting
figure(1)
plot(xArray(2,:),xArray(1,:),':','LineWidth',10);
hold on;
plot(data_log(:,4),data_log(:,5),'-k','LineWidth',10);
hold on;
plot(x_mArray(2,:),x_mArray(1,:),'-.','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',12)
grid on;
grid minor;
legend('WyNDA','Experiment','3D Model')
ylabel('y (m)','FontSize',24)
xlabel('x (m)','FontSize',24)
ylim([-1 1.5])
xlim([-0.5 2])