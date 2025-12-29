function generate_channel_param(K, M, user_center, radius, GBS_pos, LBS_pos, UAV_pos, Jammer_pos, output_file)
    % 参数说明:
    % K: 用户数量
    % M: 总路径数（包括直射路径和散射路径）
    % user_center: 用户区域中心坐标 [x, y, z]
    % radius: 用户区域的半径
    % GBS_pos: GBS 的位置 [x, y, z]
    % LBS_pos: LBS 的位置 [x, y, z]
    % UAV_pos: 无人机的位置 [x, y, z]
    % Jammer_pos: 干扰源的位置 [x, y, z]
    % output_file: 输出文件名（.dat 文件）

    % 生成 K 个用户的位置
    theta_users = 2 * pi * rand(1, K); % 随机生成角度
    r_users = radius * sqrt(rand(1, K)); % 随机生成半径（均匀分布在圆内）
    user_x = user_center(1) + r_users .* cos(theta_users); % 用户 x 坐标
    user_y = user_center(2) + r_users .* sin(theta_users); % 用户 y 坐标
    user_z = user_center(3) * ones(1, K); % 用户 z 坐标
    user_positions = [user_x; user_y; user_z]; % 用户位置矩阵 (3 x K)

    % 初始化存储变量
    % 基站到用户的角度和距离
    theta_G_t = zeros(M, K); % GBS 到用户的发射仰角（包括直射和散射）
    phi_G_t = zeros(M, K); % GBS 到用户的发射方位角
    theta_L_t = zeros(M, K); % LBS 到用户的发射仰角
    phi_L_t = zeros(M, K); % LBS 到用户的发射方位角
    d_single_G = zeros(K, 1); % GBS 到用户的距离（仅直射路径）
    d_single_L = zeros(K, 1); % LBS 到用户的距离（仅直射路径）

    % 基站到无人机的角度和距离
    theta_GA_t = zeros(M, 1); % GBS 到无人机的发射仰角
    phi_GA_t = zeros(M, 1); % GBS 到无人机的发射方位角
    theta_LA_t = zeros(M, 1); % LBS 到无人机的发射仰角
    phi_LA_t = zeros(M, 1); % LBS 到无人机的发射方位角
    d_GA = zeros(1, 1); % GBS 到无人机的距离（仅直射路径）
    d_LA = zeros(1, 1); % LBS 到无人机的距离（仅直射路径）

    % 干扰源到无人机的角度和距离
    theta_JA_t = zeros(M, 1); % 干扰源到无人机的发射仰角
    phi_JA_t = zeros(M, 1); % 干扰源到无人机的发射方位角
    d_JA = zeros(1, 1); % 干扰源到无人机的距离（仅直射路径）

    % 无人机接收端的角度信息
    theta_GA_r = zeros(M, 1); % 无人机从 GBS 接收的仰角
    phi_GA_r = zeros(M, 1); % 无人机从 GBS 接收的方位角
    theta_LA_r = zeros(M, 1); % 无人机从 LBS 接收的仰角
    phi_LA_r = zeros(M, 1); % 无人机从 LBS 接收的方位角
    theta_JA_r = zeros(M, 1); % 无人机从干扰源接收的仰角
    phi_JA_r = zeros(M, 1); % 无人机从干扰源接收的方位角

    % 计算基站到用户的角度和距离（包括直射和散射路径）
    for k = 1:K
        user_pos = user_positions(:, k);

        % GBS 到用户
        diff_G = user_pos - GBS_pos;
        theta_G_t(1, k) = atan2(diff_G(3), norm(diff_G(1:2))); % 直射路径仰角
        phi_G_t(1, k) = atan2(diff_G(2), diff_G(1)); % 直射路径方位角
        d_single_G(k) = norm(diff_G); % 直射路径距离

        % LBS 到用户
        diff_L = user_pos - LBS_pos;
        theta_L_t(1, k) = atan2(diff_L(3), norm(diff_L(1:2))); % 直射路径仰角
        phi_L_t(1, k) = atan2(diff_L(2), diff_L(1)); % 直射路径方位角
        d_single_L(k) = norm(diff_L); % 直射路径距离

        % 散射路径
        for m = 2:M
            % 随机生成散射点（在直射路径的周围）
            scatter_pos = GBS_pos + 0.5 * (user_pos - GBS_pos); % 直射路径上的点
            scatter_pos = scatter_pos + 10 * randn(3, 1); % 在路径周围随机偏移
            diff_G_scatter = scatter_pos - GBS_pos;
            theta_G_t(m, k) = atan2(diff_G_scatter(3), norm(diff_G_scatter(1:2))); % 散射路径仰角
            phi_G_t(m, k) = atan2(diff_G_scatter(2), diff_G_scatter(1)); % 散射路径方位角

            scatter_pos = LBS_pos + 0.5 * (user_pos - LBS_pos); % 直射路径上的点
            scatter_pos = scatter_pos + 10 * randn(3, 1); % 在路径周围随机偏移
            diff_L_scatter = scatter_pos - LBS_pos;
            theta_L_t(m, k) = atan2(diff_L_scatter(3), norm(diff_L_scatter(1:2))); % 散射路径仰角
            phi_L_t(m, k) = atan2(diff_L_scatter(2), diff_L_scatter(1)); % 散射路径方位角
        end
    end

    % 计算基站到无人机的角度和距离（包括直射和散射路径）
    diff_GA = UAV_pos - GBS_pos;
    theta_GA_t(1) = atan2(diff_GA(3), norm(diff_GA(1:2))); % 直射路径仰角
    phi_GA_t(1) = atan2(diff_GA(2), diff_GA(1)); % 直射路径方位角
    d_GA = norm(diff_GA); % 直射路径距离

    diff_LA = UAV_pos - LBS_pos;
    theta_LA_t(1) = atan2(diff_LA(3), norm(diff_LA(1:2))); % 直射路径仰角
    phi_LA_t(1) = atan2(diff_LA(2), diff_LA(1)); % 直射路径方位角
    d_LA = norm(diff_LA); % 直射路径距离

    % 散射路径
    for m = 2:M
        % 随机生成散射点（在直射路径的周围）
        scatter_pos = GBS_pos + 0.5 * (UAV_pos - GBS_pos); % 直射路径上的点
        scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
        diff_GA_scatter = scatter_pos - GBS_pos;
        theta_GA_t(m) = atan2(diff_GA_scatter(3), norm(diff_GA_scatter(1:2))); % 散射路径仰角
        phi_GA_t(m) = atan2(diff_GA_scatter(2), diff_GA_scatter(1)); % 散射路径方位角

        scatter_pos = LBS_pos + 0.5 * (UAV_pos - LBS_pos); % 直射路径上的点
        scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
        diff_LA_scatter = scatter_pos - LBS_pos;
        theta_LA_t(m) = atan2(diff_LA_scatter(3), norm(diff_LA_scatter(1:2))); % 散射路径仰角
        phi_LA_t(m) = atan2(diff_LA_scatter(2), diff_LA_scatter(1)); % 散射路径方位角
    end

    % 计算干扰源到无人机的角度和距离（包括直射和散射路径）
    diff_JA = UAV_pos - Jammer_pos;
    theta_JA_t(1) = atan2(diff_JA(3), norm(diff_JA(1:2))); % 直射路径仰角
    phi_JA_t(1) = atan2(diff_JA(2), diff_JA(1)); % 直射路径方位角
    d_JA = norm(diff_JA); % 直射路径距离

    % 散射路径
    for m = 2:M
        % 随机生成散射点（在直射路径的周围）
        scatter_pos = Jammer_pos + 0.5 * (UAV_pos - Jammer_pos); % 直射路径上的点
        scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
        diff_JA_scatter = scatter_pos - Jammer_pos;
        theta_JA_t(m) = atan2(diff_JA_scatter(3), norm(diff_JA_scatter(1:2))); % 散射路径仰角
        phi_JA_t(m) = atan2(diff_JA_scatter(2), diff_JA_scatter(1)); % 散射路径方位角
    end

    % 计算无人机接收端的角度信息（包括直射和散射路径）
    for m = 1:M
        % GBS 到无人机的接收角度
        if m == 1
            diff_GA_r = GBS_pos - UAV_pos; % 直射路径
        else
            scatter_pos = GBS_pos + 0.5 * (UAV_pos - GBS_pos); % 直射路径上的点
            scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
            diff_GA_r = scatter_pos - UAV_pos; % 散射路径
        end
        theta_GA_r(m) = atan2(diff_GA_r(3), norm(diff_GA_r(1:2))); % 接收仰角
        phi_GA_r(m) = atan2(diff_GA_r(2), diff_GA_r(1)); % 接收方位角

        % LBS 到无人机的接收角度
        if m == 1
            diff_LA_r = LBS_pos - UAV_pos; % 直射路径
        else
            scatter_pos = LBS_pos + 0.5 * (UAV_pos - LBS_pos); % 直射路径上的点
            scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
            diff_LA_r = scatter_pos - UAV_pos; % 散射路径
        end
        theta_LA_r(m) = atan2(diff_LA_r(3), norm(diff_LA_r(1:2))); % 接收仰角
        phi_LA_r(m) = atan2(diff_LA_r(2), diff_LA_r(1)); % 接收方位角

        % 干扰源到无人机的接收角度
        if m == 1
            diff_JA_r = Jammer_pos - UAV_pos; % 直射路径
        else
            scatter_pos = Jammer_pos + 0.5 * (UAV_pos - Jammer_pos); % 直射路径上的点
            scatter_pos = scatter_pos + radius * randn(3, 1); % 在路径周围随机偏移
            diff_JA_r = scatter_pos - UAV_pos; % 散射路径
        end
        theta_JA_r(m) = atan2(diff_JA_r(3), norm(diff_JA_r(1:2))); % 接收仰角
        phi_JA_r(m) = atan2(diff_JA_r(2), diff_JA_r(1)); % 接收方位角
    end

    save(output_file, 'theta_G_t', 'phi_G_t', 'theta_L_t', 'phi_L_t', ...
         'theta_GA_t', 'phi_GA_t', 'theta_LA_t', 'phi_LA_t', ...
         'theta_JA_t', 'phi_JA_t', 'd_single_G', 'd_single_L', ...
         'd_GA', 'd_LA', 'd_JA', 'theta_GA_r', 'phi_GA_r', ...
         'theta_LA_r', 'phi_LA_r', 'theta_JA_r', 'phi_JA_r');
end