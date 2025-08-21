% Bounds [lambda_v, lambda_t, Rx, Rt, Px, Pt]
lb = [0.7, 0.9, 0.01, 0.01, 0.01, 0.01];
ub = [0.95, 0.995, 1e5, 1e5, 1e5, 1e5];
nvars = 6;

% Parameter ground truth
m_true = 0.85613;
lf_true = 0.06874;
lr_true = 0.06726;
I_true = 0.00794;
l_true = lf_true+lr_true;

opts = optimoptions('particleswarm', ...
    'Display', 'iter', ...
    'SwarmSize', 150);


best_param = particleswarm(@(p) wynda_wmr_costfunct(p, data_log, n, r, m_true, lf_true, lr_true, I_true, l_true), nvars, lb, ub, opts)