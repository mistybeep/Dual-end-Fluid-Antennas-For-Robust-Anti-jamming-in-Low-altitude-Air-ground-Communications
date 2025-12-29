function [v_A, mu] = find_vA_rA(H_GA, H_LA, W_G, w_L, v_A, P_J, Sigma_JA, Sigma_LA, Sigma_GA, ...
                                                     G_GA, G_LA, G_JA, norm_A, t_G, t_L, ...
                                                     Q, R, data, err_range, sigma2, t_pos_J, r_A, lambda)
    % mu = ones(Q,1)/Q;sum_A, F_LA, H_GA, H_LA
    % nu = ones(R,1)/R;

    mu = ones(Q*R,1)/(Q*R);

    N_A = size(H_LA, 1);
    % N_J = JA.N;
    theta_JA_t = data.theta_JA_t;
    phi_JA_t = data.phi_JA_t;
    theta_JA_r = data.theta_JA_r;
    phi_JA_r = data.phi_JA_r;
    sum_A = zeros(N_A, N_A);
    A = zeros(N_A, N_A, Q*R);
    for idx = 1:Q*R
        [q, r] = ind2sub([Q, R], idx);
        theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
        phi_r   = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
        [H_JA, ~, ~] = generate_channel(t_pos_J, r_A, lambda, theta_JA_t, phi_JA_t, theta_q, phi_r, Sigma_JA, false);
        A(:,:,idx) = H_JA * H_JA';
        sum_A = sum_A + mu(idx)*A(:,:,idx);
    end
    % for q = 1:Q
    %     theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
    %     for r = 1:R
    %         phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
    %         [H_JA, ~, ~] = generate_channel(t_pos_J, r_pos_multi, lambda, theta_JA_t, phi_JA_t, theta_q, phi_r, Sigma_JA, false);
    %         A(:,:,(q-1)*R + r) = H_JA * H_JA';
    %     end
    % end
    v_A = inv(H_LA*(w_L*w_L')*H_LA' + H_GA*(W_G*W_G')*H_GA' + P_J*sum_A + sigma2*eye(N_A)) * H_LA*w_L;
    v_A = v_A / norm(v_A, 2);
    g = zeros(2,1);
    g(2) = -1000;
    max_iter = 1000;
    s = 0;
    while(abs(g(1) - g(2))>1e-4)
        s = s + 1;
        if s>max_iter
            break;
        end
        g(1) = g(2);

        numerator = zeros(Q*R, 1); % 存储每个 q 的分子
        denominator = 0; % 存储分母
        for i = 1:Q*R
            numerator(i) = v_A' * A(:, :, i) * v_A;
            denominator = denominator + numerator(i);
        end
        mu = numerator / denominator;

        v_A = inv(H_LA*(w_L*w_L')*H_LA' + H_GA*(W_G*W_G')*H_GA' + P_J*sum_A + sigma2*eye(N_A)) * H_LA*w_L;
        v_A = v_A / norm(v_A, 2);

        sum_A = zeros(N_A, N_A);
        for idx = 1:Q*R
            sum_A = sum_A + mu(idx)*A(:,:,idx);
        end
        v_A = inv(H_LA*(w_L*w_L')*H_LA' + H_GA*(W_G*W_G')*H_GA' + P_J*sum_A + sigma2*eye(N_A)) * H_LA*w_L;
        v_A = v_A / norm(v_A, 2);

        g(2) = abs(v_A'*H_LA*w_L)^2 / (sum_square_abs(v_A' * H_GA * W_G) + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2);

    end

end