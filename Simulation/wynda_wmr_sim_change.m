%% ============================================================
%  Research Code
%  Authors : Agus Hasan & Shinta Kosmaga
%  Purpose : WMR piecwise time-varying parameter estimation using WyNDA
%% ============================================================
clear all; close all;

%% Number of variables and coefficients
n = 5;
r = 20;

%% Time horizon
tf  = 120;
dt  = 0.05;
t   = dt:dt:tf;
nsteps = length(t);

%% System description
A = eye(n);
C = eye(n);

%% Noise
R = 0;

%% State and parameter initialization
x        = zeros(n,1);
xbar     = x;
y        = x;
thetabar = zeros(r,1);
Ibar     = zeros(1,1);
dfbar    = zeros(1,1);
drbar    = zeros(1,1);
mbar     = zeros(1,1);
ebar     = xbar-x;

% Initial parameter
mass1 = 0.85613;  df1 = 0.06874;  dr1 = 0.06726;  Kf1   = 0.04;  Kr1 = 0.0435;   I1  = 0.00794;
mass2 = 0.94218;  df2 = 0.0693;   dr2 = 0.0667;   Kf2 = 0.04;    Kr2 = 0.0435;   I2 = 0.0083;
mass3 = 1.11228;  df3 = 0.07859;  dr3 = 0.05741;  Kf3 = 0.04;    Kr3 = 0.0435;   I3 = 0.01008;
mass4 = 1.28419;  df4 = 0.08275;  dr4 = 0.05325;  Kf4 = 0.04;    Kr4 = 0.0435;   I4 = 0.01035;

%% Scenario
Scenario(1).mass = [mass1 mass2];
Scenario(1).df   = [df1   df2];
Scenario(1).dr   = [dr1   dr2];
Scenario(1).Kf   = [Kf1   Kf2];
Scenario(1).Kr   = [Kr1   Kr2];
Scenario(1).I    = [I1    I2];

Scenario(3).mass = [mass1 mass3];
Scenario(3).df   = [df1   df3];
Scenario(3).dr   = [dr1   dr3];
Scenario(3).Kf   = [Kf1   Kf3];
Scenario(3).Kr   = [Kr1   Kr3];
Scenario(3).I    = [I1    I3];

Scenario(2).mass = [mass1 mass4];
Scenario(2).df   = [df1   df4];
Scenario(2).dr   = [dr1   dr4];
Scenario(2).Kf   = [Kf1   Kf4];
Scenario(2).Kr   = [Kr1   Kr4];
Scenario(2).I    = [I1    I4];

nScenario = length(Scenario);

%% WyNDA Hyperparameter Initialization
lambdav = 0.995;
lambdat = 0.9964;
Rx = 100000*eye(n);
Rt = 3.75e+04*eye(n);

%% Scenario loop
for s = 1:nScenario

    fprintf('Running Scenario %d...\n',s)

    %% Reset states
    x = zeros(n,1);
    xbar = x;
    thetabar = zeros(r,1);

    Px = 8.93e+04*eye(n);
    Pt = 100000*eye(r);
    Gamma = zeros(n,r);

    %% initial control inputs
    S = 0;    V = 0;    B = 0;

    %% Data storage
    SArray         = []; VArray         = []; BArray         = [];
    xArray         = []; xbarArray      = []; yArray         = []; thetabarArray  = []; ebarArray      = [];
    IbarArray      = []; dfbarArray     = []; drbarArray     = []; mbarArray      = [];
    KfArray        = []; KrArray        = []; IArray         = []; dfArray        = []; drArray        = []; mArray         = [];

    %% True WMR parameters
    for i = 1:nsteps
        % parameter change at 35s
        if i < 700
            idx = 1;
        else
            idx = 2;
        end
        
        mArray(i) = Scenario(s).mass(idx);
        dfArray(i) = Scenario(s).df(idx);
        drArray(i) = Scenario(s).dr(idx);
        KfArray(i) = Scenario(s).Kf(idx);
        KrArray(i) = Scenario(s).Kr(idx);
        IArray(i) = Scenario(s).I(idx);
        
    end
    
    %% WyNDA simulation loop
    for i=1:nsteps
        % parameter change at 35s
        if i < 700
            idx = 1;
        else
            idx = 2;
        end

        % True WMR Parameter
        mArray(i) = Scenario(s).mass(idx);
        dfArray(i) = Scenario(s).df(idx);
        drArray(i) = Scenario(s).dr(idx);
        KfArray(i) = Scenario(s).Kf(idx);
        KrArray(i) = Scenario(s).Kr(idx);
        IArray(i) = Scenario(s).I(idx);
        
        m = mArray(i);
        df = dfArray(i);
        dr = drArray(i);
        Kf = KfArray(i);
        Kr = KrArray(i);
        I = IArray(i);
        d=df+dr;

        % Stored data
        VArray         = [VArray V];
        SArray         = [SArray S];
        BArray         = [BArray B];
        xArray         = [xArray x];
        xbarArray      = [xbarArray xbar];      
        yArray         = [yArray y];
        thetabarArray  = [thetabarArray thetabar]; 
        IbarArray      = [IbarArray Ibar];
        dfbarArray     = [dfbarArray dfbar];
        drbarArray     = [drbarArray drbar];
        mbarArray      = [mbarArray mbar];
        ebarArray      = [ebarArray ebar];
    
        % Input profile
        if i<200
            V = 0.1;
            S = -0.15;
        elseif i>=200 && i<800
            V = 0.05;
            S = 0.05;
        elseif i>=800 && i<1100
            V = 0.12;
            S = 0.25;
        else
            V = 0.1;
            S = 0.05;
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

    %% Final Estimation Statistics (Steady-State)
    Result(s).x_true = xArray;
    Result(s).x_est = xbarArray;
    Result(s).err = ebarArray;
    
    % True Parameter
    Result(s).m_true = mArray;
    Result(s).df_true = dfArray;
    Result(s).dr_true = drArray;
    Result(s).I_true = IArray;

    % Estimated Parameter
    Result(s).m_est = mbarArray;
    Result(s).df_est = dfbarArray;
    Result(s).dr_est = drbarArray;
    Result(s).I_est = IbarArray;
    
    % Average (Steady state)
    Result(s).mbar_avg1 = mean(mbarArray(550:700));
    Result(s).Ibar_avg1 = mean(IbarArray(550:700));
    Result(s).dfbar_avg1 = mean(dfbarArray(550:700));
    Result(s).drbar_avg1 = mean(drbarArray(550:700));
    
    Result(s).mbar_avg2 = mean(mbarArray(end-150:end));
    Result(s).Ibar_avg2 = mean(IbarArray(end-150:end));
    Result(s).dfbar_avg2 = mean(dfbarArray(end-150:end));
    Result(s).drbar_avg2 = mean(drbarArray(end-150:end));
    
    % RMSE
    Result(s).mbarerr1 = sqrt(mean((mbarArray(550:700) - Scenario(s).mass(1)).^2));
    Result(s).dfbarerr1 = sqrt(mean((dfbarArray(550:700) - Scenario(s).df(1)).^2));
    Result(s).drbarerr1 = sqrt(mean((drbarArray(550:700) - Scenario(s).dr(1)).^2));
    Result(s).Iarerr1 = sqrt(mean((IbarArray(550:700) - Scenario(s).I(1)).^2));
    
    Result(s).mbarerr2 = sqrt(mean((mbarArray(end-150:end) - Scenario(s).mass(2)).^2));
    Result(s).dfbarerr2 = sqrt(mean((dfbarArray(end-150:end) - Scenario(s).df(2)).^2));
    Result(s).drbarerr2 = sqrt(mean((drbarArray(end-150:end) - Scenario(s).dr(2)).^2));
    Result(s).Iarerr2 = sqrt(mean((IbarArray(end-150:end) - Scenario(s).I(2)).^2));

    % MAPE
    Result(s).mbarerrp1 = mean(abs((Scenario(s).mass(1)-mbarArray(550:700))./Scenario(s).mass(1))) * 100;
    Result(s).dfbarerrp1 = mean(abs((Scenario(s).df(1)-dfbarArray(550:700))./Scenario(s).df(1))) * 100;
    Result(s).drbarerrp1 = mean(abs((Scenario(s).dr(1)-drbarArray(550:700))./Scenario(s).dr(1))) * 100;
    Result(s).Iarerrp1 = mean(abs((Scenario(s).I(1)-IbarArray(550:700))./Scenario(s).I(1))) * 100;

    Result(s).mbarerrp2 = mean(abs((Scenario(s).mass(2)-mbarArray(end-150:end))./Scenario(s).mass(2))) * 100;
    Result(s).dfbarerrp2 = mean(abs((Scenario(s).df(2)-dfbarArray(end-150:end))./Scenario(s).df(2))) * 100;
    Result(s).drbarerrp2 = mean(abs((Scenario(s).dr(2)-drbarArray(end-150:end))./Scenario(s).dr(2))) * 100;
    Result(s).Iarerrp2 = mean(abs((Scenario(s).I(2)-IbarArray(end-150:end))./Scenario(s).I(2))) * 100;
end

%% Plotting

% Input
figure(1)
subplot(2,1,1)
plot(t, SArray, '-k', 'LineWidth', 10)
formatAxis('\gamma (rad)', 't (s)')

subplot(2,1,2)
plot(t, VArray, '-k', 'LineWidth', 10)
formatAxis('v (m/s)', 't (s)')


% Parameter estimation
figure(2)
colors = lines(nScenario);
paramName = {'m','I','df','dr'};
yLabel    = {'m (kg)','Iy (kgm^2)','df (m)','dr (m)'};
yLimit    = {[0 2],[0 0.02],[0 0.1],[0 0.1]};

for p = 1:4
    subplot(2,2,p); hold on;
    
    for s = 1:nScenario
        plot(t, Result(s).([paramName{p} '_true']), '-', ...
            'Color', colors(s,:), 'LineWidth', 5, ...
            'DisplayName', ['Scenario ' num2str(s)]);
        
        plot(t, Result(s).([paramName{p} '_est']), ':', ...
            'Color', colors(s,:), 'LineWidth', 5, ...
            'HandleVisibility','off');
    end
    
    formatAxis(yLabel{p}, 't (s)')
    ylim(yLimit{p})
    
    if p == 1
        legend('Location','best')
    end
end

function formatAxis(ylabelText, xlabelText)
    set(gca,'Color','white','LineWidth',2,'FontSize',12)
    grid on; grid minor;
    ylabel(ylabelText,'FontSize',18)
    xlabel(xlabelText,'FontSize',18)
end
