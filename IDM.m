function ahat = IDM(vi, delta_v, delta_d)
    vf= 30;
    A = 3;
    b = 3.2735;
    s0= 2;
    T = 2.2;
    s_star = s0 + max(0, (vi*T + (vi * delta_v) / (2 * sqrt(A*b))) );
    epsilon = 1e-20;
    ahat = A*(1 - (vi/vf)^4 - (s_star/(delta_d + epsilon))^2);
end