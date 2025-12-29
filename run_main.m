% 主脚本

num_simulations = 1; % 仿真次数
rates = zeros(num_simulations, 11); % 存储每次仿真的结果

lambda = 0.01;
norm_A = 4;
M = 8;
N_G = 12;
N_L = 12;
N_J = 12;
N_A = 4;
K = 4;
P_dBm = 20;
P_J_dBm = 30;
sigma2_dBm = -80;
gamma = 1 * ones(K, 1);
err_range = 4 / 180 * pi;
alpha = 2.8;
rho_0_dB = -40;
Q = 6;
R = 6;
data = load_dataset('dataset.mat');

% 运行仿真
for sim = 1:num_simulations %20
    fprintf('运行仿真 %d/%d...\n', sim, num_simulations);
    result = run_simulation(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R);
    if ~any(isnan(result)) % 检查是否包含 NaN
        rates(sim, :) = result; % 如果不包含 NaN，存入 rates
    else
        rates(sim, :) = []; % 如果包含 NaN，跳过该结果
    end
end

% 计算平均速率
average_rate = mean(rates, 1);
figure;
plot(0:10,average_rate,'r-o')
ylim([0 max(average_rate)+0.1]);
fprintf('平均速率: %.4f\n', average_rate);

% load("RATE.mat");
% load("MA_RATE.mat");
% plot(1:20,RATE)
% hold on;
% plot(1:20,MA_RATE);
% ylim([0 2.5]);

function rate = run_simulation(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R)
    rng(42)
    P = 10^(P_dBm / 10) * 1e-3; % 转换为 W
    P_J = 10^(P_J_dBm / 10) * 1e-3; % 转换为 W
    sigma2 = 10^(sigma2_dBm / 10) * 1e-3; % 转换为 W
    rho_0 = 10^(rho_0_dB / 10); % 转换为线性单位
    
    % user_center = [40; 30; 0]; % 用户区域中心
    % radius = 10; % 用户区域半径
    % GBS_pos = [0; 0; 30]; % GBS 位置
    % LBS_pos = [0; 100; 30]; % LBS 位置
    % UAV_pos = [10; 60; 100]; % 无人机位置
    % Jammer_pos = [40; 60; 0]; % 干扰源位置
    
    % t_G = [0:lambda/2:(N_G-1)*lambda/2; zeros(1, N_G)];
    factors = [];
    for i = 1:sqrt(N_G)
        if mod(N_G, i) == 0
            factors = [factors; i, N_G/i];
        end
    end

    % 选择 N_L2 最小的分解（即 N_L1 最大的情况）
    if ~isempty(factors)
        [N_G1, N_G2] = deal(factors(end, 1), factors(end, 2));
    else
        % 如果 N_L 是质数，退化成 ULA（线性阵列）
        N_G1 = 1;
        N_G2 = N_G;
    end
    x_pos = linspace(-(N_G1-1)*lambda/4, (N_G1-1)*lambda/4, N_G1); % x 方向
    y_pos = linspace(-(N_G2-1)*lambda/4, (N_G2-1)*lambda/4, N_G2); % y 方向
    % 生成网格坐标
    [X, Y] = meshgrid(x_pos, y_pos);
    % 组合成 2×N_L 的矩阵 [x; y]
    t_G = [X(:)'; Y(:)'];


    factors = [];
    for i = 1:sqrt(N_L)
        if mod(N_L, i) == 0
            factors = [factors; i, N_L/i];
        end
    end

    % 选择 N_L2 最小的分解（即 N_L1 最大的情况）
    if ~isempty(factors)
        [N_L1, N_L2] = deal(factors(end, 1), factors(end, 2));
    else
        % 如果 N_L 是质数，退化成 ULA（线性阵列）
        N_L1 = 1;
        N_L2 = N_L;
    end
    x_pos = linspace(-(N_L1-1)*lambda/4, (N_L1-1)*lambda/4, N_L1); % x 方向
    y_pos = linspace(-(N_L2-1)*lambda/4, (N_L2-1)*lambda/4, N_L2); % y 方向
    % 生成网格坐标
    [X, Y] = meshgrid(x_pos, y_pos);
    % 组合成 2×N_L 的矩阵 [x; y]
    t_L = [X(:)'; Y(:)'];
    % t_L = generate_upa_coordinates(norm_A, N_L, lambda);
    % load t_L.mat t_L;
    t_J = [0:lambda/2:(N_J-1)*lambda/2; zeros(1, N_J)];
    r_pos_single = zeros(2, K); % 单天线用户位置（假设在原点）
    % r_A = generate_upa_coordinates(norm_A, N_A, lambda);
    load r_A.mat r_A;
    % factors = [];
    % for i = 1:sqrt(N_A)
    %     if mod(N_A, i) == 0
    %         factors = [factors; i, N_A/i];
    %     end
    % end
    % 
    % % 选择 N_L2 最小的分解（即 N_L1 最大的情况）
    % if ~isempty(factors)
    %     [N_A1, N_A2] = deal(factors(end, 1), factors(end, 2));
    % else
    %     % 如果 N_L 是质数，退化成 ULA（线性阵列）
    %     N_A1 = 1;
    %     N_A2 = N_A;
    % end
    % x_pos = linspace(-(N_A1-1)*lambda/4, (N_A1-1)*lambda/4, N_A1); % x 方向
    % y_pos = linspace(-(N_A2-1)*lambda/4, (N_A2-1)*lambda/4, N_A2); % y 方向
    % % 生成网格坐标
    % [X, Y] = meshgrid(x_pos, y_pos);
    % % 组合成 2×N_L 的矩阵 [x; y]
    % r_A = [X(:)'; Y(:)'];
    
    % load('./dataset.mat');
    % Sigma_GA_Test = zeros(100,M,M);
    % Sigma_LA_Test = zeros(100,M,M);
    % Sigma_JA_Test = zeros(100,M,M);
    % Sigma_GK_Test = zeros(100,M,M,K);
    % Sigma_LK_Test = zeros(100,M,M,K);
    load Sigma.mat Sigma_GA_Test Sigma_GK_Test Sigma_JA_Test Sigma_LA_Test Sigma_LK_Test
    % load("./Sigma.mat");
    Sigma_GA = squeeze(Sigma_GA_Test(sim,:,:));
    Sigma_LA = squeeze(Sigma_LA_Test(sim,:,:));
    Sigma_JA = squeeze(Sigma_JA_Test(sim,:,:));
    Sigma_GK = squeeze(Sigma_GK_Test(sim, :, :, :));
    Sigma_LK = squeeze(Sigma_LK_Test(sim, :, :, :));
    % Sigma_GA = generate_Sigma(M, rho_0, alpha, data.d_GA);
    % Sigma_LA = generate_Sigma(M, rho_0, alpha, data.d_LA);
    % Sigma_JA = generate_Sigma(M, rho_0, alpha, data.d_JA);
    % Sigma_GA_Test = zeros(200,M,M);
    % Sigma_LA_Test = zeros(200,M,M);
    % Sigma_JA_Test = zeros(200,M,M);
    % Sigma_GK_Test = zeros(200,M,M,K);
    % Sigma_LK_Test = zeros(200,M,M,K);
    % for i = 1:200
    %     Sigma_GA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_GA);
    %     Sigma_LA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_LA);
    %     Sigma_JA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_JA);
    %     for k = 1:K
    %         Sigma_GK_Test(i,:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_G(k));
    %         Sigma_LK_Test(i,:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_L(k));
    %     end
    % end
    % save("Sigma.mat", 'Sigma_GA_Test', 'Sigma_LA_Test', 'Sigma_JA_Test', 'Sigma_GK_Test','Sigma_LK_Test')
    % 生成基站 G 的信道
    H_G = zeros(K, N_G); % 基站 G 到 K 个单天线用户的信道 (K x Nt)
    % Sigma_GK = zeros(M, M, K);
    % Sigma_GK = squeeze(Sigma_GK_Test(sim, :, :, :));
    for k = 1:K
        r_pos_k = r_pos_single(:, k); % 第 k 个单天线用户的位置 (2 x 1)
        % Sigma_GK(:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_G(k));
        H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, data.theta_G_t(:, k), data.phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
    end
    [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
    % 生成基站 L 的信道
    H_L = zeros(K, N_L); % 基站 L 到 K 个单天线用户的信道 (K x Nt)
    % Sigma_LK = zeros(M, M, K);
    % Sigma_LK = squeeze(Sigma_LK_Test(sim, :, :, :));
    for k = 1:K
        r_pos_k = r_pos_single(:, k); % 第 k 个单天线用户的位置 (2 x 1)
        % Sigma_LK(:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_L(k));
        H_L(k, :) = generate_channel(t_L, r_pos_k, lambda, data.theta_L_t(:, k), data.phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
    end
    % H_LA = generate_channel(t_L, r_A, lambda, theta_LA_t, phi_LA_t, theta_LA_r, phi_LA_r, Sigma_LA, false);
    % [H_GA, ~, G_GA] = generate_channel(t_G, r_A, lambda, theta_GA_t, phi_GA_t, theta_GA_r, phi_GA_r, Sigma_GA, false);
    [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
    
    [~, ~, G_JA] = generate_channel(t_J, r_A, lambda, data.theta_JA_t, data.phi_JA_t, data.theta_JA_r, data.phi_JA_r, Sigma_JA, false);
    
    w_L = randn(N_L, 1) + 1j * randn(N_L, 1);
    w_L = w_L * sqrt(P) / norm(w_L, 2);
    v_A = randn(N_A, 1) + 1j*randn(N_A, 1);
    v_A = v_A / norm(v_A, 2);
    
    W_G = initialize_W_G(H_G, H_L, w_L, gamma, sigma2, K);
    
    L_max = 10;

    mu = ones(Q,1)/Q;
    nu = ones(R,1)/R;
    rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range);
    r_A_org = r_A;
    t_L_org = t_L;
    rate = rate_new;
    % plot_channel_gain(data, t_J, r_A, Sigma_JA, lambda, norm_A)

    for l = 1:L_max
        y = find_y(H_G, H_L, W_G, w_L, sigma2, K);
    
        z = find_z(H_G, H_L, W_G, w_L, sigma2, y, K);
        [v_A, mu, nu, ~] = find_v_A(H_GA, H_LA, W_G, w_L, v_A, P_J, Sigma_JA, mu, nu, Q, R, data, err_range, sigma2, t_J, r_A, lambda);

        w_L = find_w_L_2(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2);
        
        W_G = find_W_G_3(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);
        
        r_A = update_r_A_PPO(r_A, data.theta_GA_r, data.phi_GA_r, Sigma_GA, ...
                        data.theta_LA_r, data.phi_LA_r, Sigma_LA, ...
                        data.theta_JA_r, data.phi_JA_r, Sigma_JA, P_J, ...
                        v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                        mu, nu, lambda, err_range, norm_A, sigma2);

        [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
        [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
        % % % % 
        % t_L = update_t_L(t_L, data.theta_LA_t, data.phi_LA_t, data.theta_G_t, data.phi_G_t, ...
        %                 Sigma_LA, Sigma_LK, F_LA, v_A, w_L, ...
        %                 y, z, gamma, norm_A, lambda, sigma2);
        % [H_LA, ~, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
        % for k = 1:K
        %     r_pos_k = r_pos_single(:, k);
        %     [H_L(k, :), ~, ~] = generate_channel(t_L, r_pos_k, lambda, data.theta_L_t(:, k), data.phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
        %     % H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, theta_G_t(:, k), phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
        % end

        rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range);
        rate = [rate rate_new];
    end
    % plot_Rec_HeatMap(v_A, r_A, N_A, lambda);
    plot_channel_gain(data, t_J, r_A, r_A_org, Sigma_JA, lambda, norm_A);
    % plot_channel_gain_user(data, t_L, r_A, r_A_org, Sigma_LA, lambda, norm_A)
    % plot_channel_gain_Tx(data, t_L, r_A, t_L_org, Sigma_LA, lambda, norm_A)
    % plot_UPA_heatmap(W_G, t_L, N_L, lambda);
end