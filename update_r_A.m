function r_A = update_r_A(r_A, theta_GA, phi_GA, Sigma_GA, ...
                            theta_LA, phi_LA, Sigma_LA, ...
                            theta_JA, phi_JA, Sigma_JA, P_J, ...
                            v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                            mu, nu, lambda, err_range, normA, sigma2)

    N_A = size(v_A, 1);
    M_GA = size(theta_GA, 1);
    M_LA = size(theta_LA, 1);
    M_JA = size(theta_JA, 1);

    Q = size(mu, 1);
    R = size(nu, 1);
    t = 0; % 初始化迭代次数
    max_iter = 100; % 最大迭代次数
    tol = 1e-3; % 容差
    % max_oscillation_count = 100; % 允许的最大交换次数
    
    % oscillation_count = 0;
    kappa = 0;
    kappa2 = -100;

    while(abs(kappa - kappa2)>tol && t<max_iter)
        t = t + 1;

        kappa2 = kappa;
        F_GA = create_F(r_A, theta_GA, phi_GA, M_GA, N_A, lambda);
        F_LA = create_F(r_A, theta_LA, phi_LA, M_LA, N_A, lambda);
        sum_A = zeros(N_A, N_A);
        for q = 1:Q
            theta_q = theta_JA - err_range/2 + (q-1)*err_range/(Q-1);
            for r = 1:R
                phi_r = phi_JA - err_range/2 + (r-1)*err_range/(R-1);
                F_JA = create_F(r_A, theta_q, phi_r, M_JA, N_A, lambda);
                sum_A = sum_A + mu(q)*nu(r)*F_JA'*Sigma_JA*(G_JA*G_JA')*Sigma_JA'*F_JA;
            end
        end
        kappa = real(abs(v_A'*F_LA'*Sigma_LA*G_LA*w_L)^2 / ...
                     real(sum(abs(v_A' * F_GA'* Sigma_GA * G_GA * W_G).^2) ...
                     + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2));
        % disp(kappa)
        % rate = log2(1+kappa);
        % kappa = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, mu, nu, Q, R, data, sigma2, lambda, err_range);

        for n = 1:N_A
            F_GA = create_F(r_A, theta_GA, phi_GA, M_GA, N_A, lambda);
            F_LA = create_F(r_A, theta_LA, phi_LA, M_LA, N_A, lambda);
            

            [gradient, bar_f_A, bar_f_JA] = cal_gradient_f3(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                                            theta_LA, phi_LA, Sigma_LA, F_GA, F_LA, ...
                                                            theta_JA, phi_JA, Sigma_JA, P_J, ...
                                                            v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                                                            mu, nu, kappa, n, lambda, err_range);
            [~, tau] = cal_hessian_f3(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                    theta_LA, phi_LA, Sigma_LA, ...
                                    theta_JA, phi_JA, Sigma_JA, P_J, ...
                                    v_A, W_G, w_L, G_GA, G_LA, G_JA, bar_f_A, bar_f_JA, ...
                                    mu, nu, kappa, n, lambda, err_range);
            % tau = 1e-3;
            
            r_Astar = r_A(:, n) - real(gradient) / real(tau);
            if n==1
                r_A_opt = r_A(:, 2:end);
                r_An = find_pos(r_A_opt.', r_Astar.', lambda, normA, norm(real(gradient) / real(tau), 2));
                if isnan(r_An(1,1)) || isnan(r_An(1,2))
                    r_A(:, n) = r_A(:, n);
                else
                    r_A(:, n) = r_An.';
                end
            else
                r_A_opt = r_A(:, 1:n-1);
                if n<N_A
                    r_A_opt = [r_A_opt, r_A(:, n+1:end)];
                end
                r_An = find_pos(r_A_opt.', r_Astar.', lambda, normA, norm(real(gradient) / real(tau), 2));
                if isnan(r_An(1,1)) || isnan(r_An(1,2))
                    r_A(:, n) = r_A(:, n);
                else
                    r_A(:, n) = r_An.';
                end
            end

        end
    end

end


function f3 = func_3(n, r, r_A, theta_GA, phi_GA, Sigma_GA, ...
                    theta_LA, phi_LA, Sigma_LA, ...
                    theta_JA, phi_JA, Sigma_JA, P_J, ...
                    v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                    mu, nu, lambda, err_range, kappa, sigma2)
    r_A(:, n) = r;
    N_A = size(v_A, 1);
    M_GA = size(theta_GA, 1);
    M_LA = size(theta_LA, 1);
    M_JA = size(theta_JA, 1);
    F_GA = create_F(r_A, theta_GA, phi_GA, M_GA, N_A, lambda);
    F_LA = create_F(r_A, theta_LA, phi_LA, M_LA, N_A, lambda);
    Q = size(mu, 1);
    R = size(nu, 1);
    sum_A = zeros(N_A, N_A);
    for q = 1:Q
        theta_q = theta_JA - err_range/2 + (q-1)*err_range/(Q-1);
        for r = 1:R
            phi_r = phi_JA - err_range/2 + (r-1)*err_range/(R-1);
            F_JA = create_F(r_A, theta_q, phi_r, M_JA, N_A, lambda);
            sum_A = sum_A + mu(q)*nu(r)*F_JA'*Sigma_JA*(G_JA*G_JA')*Sigma_JA'*F_JA;
        end
    end
    f3 = kappa*real(sum(abs(v_A' * F_GA'* Sigma_GA * G_GA * W_G).^2) ...
        + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2) - real(abs(v_A'*F_LA'*Sigma_LA*G_LA*w_L)^2);
end