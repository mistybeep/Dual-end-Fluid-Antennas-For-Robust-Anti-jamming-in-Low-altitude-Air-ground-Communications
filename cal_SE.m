function SE = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range)
    N_A = size(v_A, 1);
    M_JA = size(Sigma_JA, 1);
    M_GA = size(Sigma_GA, 1);
    M_LA = size(Sigma_LA, 1);

    % 初始化 sum_A
    sum_A = zeros(N_A, N_A, 'like', 1+1i);  % 复数零矩阵

    % 遍历 q 和 r
    for q = 1:Q
        theta_q = data.theta_JA_r - err_range / 2 + (q-1) * err_range / (Q-1);
        for r = 1:R
            phi_r = data.phi_JA_r - err_range / 2 + (r-1) * err_range / (R-1);
            F_JA = create_F(r_A, theta_q, phi_r, M_JA, N_A, lambda);  % 调用 create_F 函数
            sum_A = sum_A + mu(q) * nu(r) * (F_JA' * Sigma_JA * (G_JA * G_JA') * Sigma_JA' * F_JA);
        end
    end

    % 计算 F_GA 和 F_LA
    F_GA = create_F(r_A, data.theta_GA_r, data.phi_GA_r, M_GA, N_A, lambda);
    F_LA = create_F(r_A, data.theta_LA_r, data.phi_LA_r, M_LA, N_A, lambda);

    % 计算分子
    numerator = abs(v_A' * F_LA' * Sigma_LA * G_LA * w_L)^2;

    % 计算分母的三项
    term1 = sum(abs(v_A' * F_GA' * Sigma_GA * G_GA * W_G).^2);
    term2 = P_J * (v_A' * sum_A * v_A);
    term3 = norm(v_A, 2)^2 * sigma2;
    denominator = term1 + term2 + term3;

    % 计算 SINR
    sinr = real(numerator / denominator);  % 取实部

    % 计算频谱效率 (SE)
    SE = log2(1 + sinr);
end