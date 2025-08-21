function cost = wynda_wmr_costfunct(param, data_log, n, r, m_true, lf_true, lr_true, I_true, l_true)
    % get parameter
    lambdav  = param(1);
    lambdat  = param(2);
    Rx_val   = param(3);
    Rt_val   = param(4);
    Px_val   = param(5);
    Pt_val   = param(6);

    % initialization
    Px = Px_val * eye(n);
    Pt = Pt_val * eye(r);
    Rx = Rx_val * eye(n);
    Rt = Rt_val * eye(n);
    Gamma = zeros(n, r);

    xbar = zeros(n,1);
    thetabar = zeros(r,1);

    mbarArray = [];
    lfbarArray = [];
    lrbarArray = [];
    IbarArray = [];

    for i=1:1:size(data_log,1)
        % get data from log
        x = zeros(n,1);
        x(1) = data_log(i,5);  % y
        x(2) = data_log(i,4);  % x
        x(3) = data_log(i,6);  % theta
        x(4) = data_log(i,7);  % y dot
        x(5) = data_log(i,8);  % theta dot
        B    = -data_log(i,10);
        V    = data_log(i,2);
        S    = -data_log(i,3);
        dt   = 0.1;
        C    = eye(n);
        A    = eye(n);
        y    = C*x;

        Phi = [V*sin(B+y(3)) 0 0 0 zeros(16,1)';
               zeros(4,1)' V*cos(B+y(3)) 0 0 0 zeros(12,1)';
               zeros(8,1)' y(5) 0 0 0 zeros(8,1)';
               zeros(12,1)' y(3) y(4)/V y(5)/V S zeros(4,1)';
               zeros(16,1)' y(3) y(4)/V y(5)/V S];

        % estimation using adaptive observer
        Kx = Px*C' / (C*Px*C' + Rx);
        Kt = Pt*Gamma'*C' / (C*Gamma*Pt*Gamma'*C' + Rt);
        Gamma = (eye(n)-Kx*C)*Gamma;

        xbar = xbar + (Kx + Gamma*Kt)*(y - C*xbar);
        thetabar = thetabar - Kt*(y - C*xbar);

        xbar = A*xbar + Phi*thetabar;
        Px = (1/lambdav)*(eye(n) - Kx*C)*Px;
        Pt = (1/lambdat)*(eye(r) - Kt*C*Gamma)*Pt;
        Gamma = Gamma - Phi;

        % estimation result
        Kfbar = 0.04;
        Krbar = 0.0435;
        mbar  = abs((2*Kfbar*dt)/thetabar(16));
        lfbar = abs((2*Krbar*l_true-(thetabar(15)*m_true/dt))/(2*Kfbar+2*Krbar));
        lrbar = l_true-lfbar;
        Ibar  = abs((2*Kfbar*lfbar*dt)/thetabar(20));

        mbarArray(end+1)  = mbar;
        lfbarArray(end+1) = lfbar;
        lrbarArray(end+1) = lrbar;
        IbarArray(end+1)  = Ibar;
    end

    % cost function
    idx_eval = size(data_log,1)-150:size(data_log);
    mbarerr  = 10*sqrt(mean((mbarArray(idx_eval) - m_true).^2));
    lfbarerr = sqrt(mean((lfbarArray(idx_eval) - lf_true).^2));
    lrbarerr = sqrt(mean((lrbarArray(idx_eval) - lr_true).^2));
    Ibarerr  = 200*sqrt(mean((IbarArray(idx_eval) - I_true).^2)) ;

    cost = mbarerr + lfbarerr + lrbarerr + Ibarerr;  % total error
end
