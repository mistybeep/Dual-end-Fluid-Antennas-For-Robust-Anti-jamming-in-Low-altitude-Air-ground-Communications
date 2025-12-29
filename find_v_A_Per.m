function [v_A, mu, nu, DeltaPRM, r_A] = find_v_A_Per(H_GA, H_LA, W_G, w_L, v_A, P_J, Sigma_LA, G_LA, t_L, Sigma_GA, G_GA, t_G, Sigma_JA, G_JA, mu, nu, Q, R, data, err_range, sigma2, t_pos_J, r_A, DeltaPRM, EpsilonPRM, lambda, norm_A)
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
        
        sum_A = zeros(N_A, N_A);
        for q = 1:Q
            theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
            for r = 1:R
                phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
                F_JA = create_F(r_A, theta_q, phi_r, size(theta_q), N_A, lambda);  % 调用 create_F 函数
                sum_A = sum_A + mu(q) * nu(r) * (F_JA' * Sigma_JA * (G_JA * G_JA') * Sigma_JA' * F_JA);
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
        % [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
        % [H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);

        sum_A = zeros(N_A, N_A);
        for q = 1:Q
            theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
            for r = 1:R
                phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
                F_JA = create_F(r_A, theta_q, phi_r, size(theta_q), N_A, lambda);  % 调用 create_F 函数
                sum_A = sum_A + mu(q) * nu(r) * (F_JA' * Sigma_JA * (G_JA * G_JA') * Sigma_JA' * F_JA);
            end
        end

        g(2) = abs(v_A'*H_LA*w_L)^2 / (sum_square_abs(v_A' * H_GA * W_G) + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2);

    end

end