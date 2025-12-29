function rate = DualEndMA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa)
    
    
    P = 10^(P_dBm / 10) * 1e-3; % 转换为 W
    P_J = 10^(P_J_dBm / 10) * 1e-3; % 转换为 W
    sigma2 = 10^(sigma2_dBm / 10) * 1e-3; % 转换为 W
    rho_0 = 10^(rho_0_dB / 10); % 转换为线性单位
    
    t_G = Gen_UPA(N_G, lambda);
    load t_L.mat t_L;
    % t_L = Gen_UPA(N_L, lambda);
    t_J = Gen_UPA(N_J, lambda);
    r_pos_single = zeros(2, K); % 单天线用户位置（假设在原点）

    load r_A.mat r_A;
    % r_A = Gen_UPA(N_A, lambda);

    load Sigma.mat Sigma_GA_Test Sigma_GK_Test Sigma_JA_Test Sigma_LA_Test Sigma_LK_Test
    Sigma_GA = squeeze(Sigma_GA_Test(sim,:,:)) / sqrt(sigma2);
    Sigma_LA = squeeze(Sigma_LA_Test(sim,:,:)) / sqrt(sigma2);
    Sigma_JA = squeeze(Sigma_JA_Test(sim,:,:)) / sqrt(sigma2);
    Sigma_GK = squeeze(Sigma_GK_Test(sim, :, :, :)) / sqrt(sigma2);
    Sigma_LK = squeeze(Sigma_LK_Test(sim, :, :, :)) / sqrt(sigma2);
    sigma2 = 1;
    % 生成基站 G 的信道
    H_G = zeros(K, N_G);
    for k = 1:K
        r_pos_k = r_pos_single(:, k);
        H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, data.theta_G_t(:, k), data.phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
    end
    [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
    
    % 生成基站 L 的信道
    H_L = zeros(K, N_L);
    for k = 1:K
        r_pos_k = r_pos_single(:, k);
        H_L(k, :) = generate_channel(t_L, r_pos_k, lambda, data.theta_L_t(:, k), data.phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
    end
    [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
    
    [~, ~, G_JA] = generate_channel(t_J, r_A, lambda, data.theta_JA_t, data.phi_JA_t, data.theta_JA_r, data.phi_JA_r, Sigma_JA, false);
    
    W_G = initialize_W_G(H_G, H_L, w_L, gamma, sigma2, K);
    
    L_max = 10;

    mu = ones(Q,1)/Q;
    nu = ones(R,1)/R;
    rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range);
    rate = rate_new;

    EpsilonPRM = norm(Sigma_JA, 'fro') * EpsKappa;

    DeltaPRM = zeros(M, M);
    
  
    for l = 1:L_max
        
        y = find_y(H_G, H_L, W_G, w_L, sigma2, K);
    
        z = find_z(H_G, H_L, W_G, w_L, sigma2, y, K);

        [v_A, mu, nu, DeltaPRM, r_A] = find_v_A(H_GA, H_LA, W_G, w_L, v_A, P_J, Sigma_LA, G_LA, t_L, Sigma_GA, G_GA, t_G, Sigma_JA, G_JA, mu, nu, Q, R, data, err_range, sigma2, t_J, r_A, DeltaPRM, EpsilonPRM, lambda, norm_A);
        [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
        [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);

        w_L = find_w_L_2(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2);
        
        W_G = find_W_G_3(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);

        t_L = update_t_L(t_L, data.theta_LA_t, data.phi_LA_t, data.theta_G_t, data.phi_G_t, ...
                        Sigma_LA, Sigma_LK, F_LA, v_A, w_L, ...
                        y, z, gamma, norm_A, lambda, sigma2);
        [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
        for k = 1:K
            r_pos_k = r_pos_single(:, k);
            [H_L(k, :), ~, ~] = generate_channel(t_L, r_pos_k, lambda, data.theta_L_t(:, k), data.phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
        end
        
        rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA+DeltaPRM, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range);
        rate = [rate rate_new];
    end
end