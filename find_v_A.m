function [v_A, mu, nu, DeltaPRM, r_A] = find_v_A(H_GA, H_LA, W_G, w_L, v_A, P_J, Sigma_LA, G_LA, t_L, Sigma_GA, G_GA, t_G, Sigma_JA, G_JA, mu, nu, Q, R, data, err_range, sigma2, t_pos_J, r_A, DeltaPRM, EpsilonPRM, lambda, norm_A)
    % mu = ones(Q,1)/Q;
    % nu = ones(R,1)/R;
    N_A = size(H_LA, 1);
    % N_J = JA.N;
    theta_JA_t = data.theta_JA_t;
    phi_JA_t = data.phi_JA_t;
    theta_JA_r = data.theta_JA_r;
    phi_JA_r = data.phi_JA_r;
    A = zeros(N_A, N_A, Q, R);
    for q = 1:Q
        theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
        for r = 1:R
            phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
            [H_JA, ~, ~] = generate_channel(t_pos_J, r_A, lambda, theta_JA_t, phi_JA_t, theta_q, phi_r, Sigma_JA, false);
            A(:,:,q,r) = H_JA * H_JA';
        end
    end

    g = zeros(2,1);
    g(2) = -100;
    max_iter = 10;
    s = 0;
    while(1)
        s = s + 1;
        if s>max_iter && g(2) > g(1)
            break;
        end
        g(1) = g(2);

        numerator = zeros(Q, 1); % 存储每个 q 的分子
        denominator = 0; % 存储分母
        for q = 1:Q
            % 计算分子部分
            sum_r = 0;
            for r = 1:R
                sum_r = sum_r + nu(r) * (v_A' * A(:, :, q, r) * v_A);
            end
            numerator(q) = sum_r;

            % 累加分母部分
            denominator = denominator + sum_r;
        end
        mu = numerator/ denominator;

        numerator = zeros(R, 1); % 存储每个 q 的分子
        denominator = 0; % 存储分母
        for r = 1:R
            % 计算分子部分
            sum_q = 0;
            for q = 1:Q
                sum_q = sum_q + mu(q) * (v_A' * A(:, :, q, r) * v_A);
            end
            numerator(r) = sum_q;

            % 累加分母部分
            denominator = denominator + sum_q;
        end
        nu = numerator/ denominator;

        DeltaPRM = UpdatePRMError(data, v_A, r_A, Sigma_JA, G_JA, mu, nu, lambda, err_range, DeltaPRM, EpsilonPRM);
        
        sum_A = zeros(N_A, N_A);
        for q = 1:Q
            theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
            for r = 1:R
                phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
                F_JA = create_F(r_A, theta_q, phi_r, size(theta_q), N_A, lambda);  % 调用 create_F 函数
                sum_A = sum_A + mu(q) * nu(r) * (F_JA' * (Sigma_JA + DeltaPRM) * (G_JA * G_JA') * (Sigma_JA + DeltaPRM)' * F_JA);
            end
        end
        v_A = (H_LA*(w_L*w_L')*H_LA' + H_GA*(W_G*W_G')*H_GA' + P_J*sum_A + sigma2*eye(N_A)) \ (H_LA*w_L);
        v_A = v_A / norm(v_A, 2);

        r_A = update_r_A(r_A, data.theta_GA_r, data.phi_GA_r, Sigma_GA, ...
                        data.theta_LA_r, data.phi_LA_r, Sigma_LA, ...
                        data.theta_JA_r, data.phi_JA_r, Sigma_JA + DeltaPRM, P_J, ...
                        v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                        mu, nu, lambda, err_range, norm_A, sigma2);
        %
        F_GA = create_F(r_A, data.theta_GA_r, data.phi_GA_r, size(data.theta_GA_r), N_A, lambda);  % 调用 create_F 函数
        H_GA = F_GA' * Sigma_GA * G_GA;
        F_LA = create_F(r_A, data.theta_LA_r, data.phi_LA_r, size(data.theta_LA_r), N_A, lambda);  % 调用 create_F 函数
        H_LA = F_LA' * Sigma_LA * G_LA;

        sum_A = zeros(N_A, N_A);
        for q = 1:Q
            theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
            for r = 1:R
                phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
                F_JA = create_F(r_A, theta_q, phi_r, size(theta_q), N_A, lambda);  % 调用 create_F 函数
                sum_A = sum_A + mu(q) * nu(r) * (F_JA' * (Sigma_JA + DeltaPRM) * (G_JA * G_JA') * (Sigma_JA + DeltaPRM)' * F_JA);
            end
        end

        g(2) = abs(v_A'*H_LA*w_L)^2 / (sum_square_abs(v_A' * H_GA * W_G) + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2);

    end

end