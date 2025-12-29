function t_L = update_t_L(t_L, theta_LA, phi_LA, theta_LK, phi_LK, ...
                            Sigma_LA, Sigma_LK, F_LA, v_A, w_L, ...
                            y, z, gamma, normA, lambda, sigma2)

    N_L = size(w_L, 1);
    M_Lk = size(Sigma_LK, 1);
    K = size(y, 1);
    M_LA = size(Sigma_LA, 1);
    t = 0;
    max_iter = 1000;
    L1 = 0;
    L2 = -100;
    tol = 1e-4;
    while(abs(L1 - L2)>tol && t<max_iter)
        t = t+1;
        L2 = L1;
        G_LA = create_G(t_L, theta_LA, phi_LA, M_LA, N_L, lambda);
        L1 = abs(v_A' * F_LA' * Sigma_LA * G_LA * w_L).^2;
        for n = 1:N_L
            % G_GK = zeros(M_Gk, N_L, K);
            grad = zeros(2, K);
            tau = zeros(K, 1);
            t_L_k_star = zeros(2, K);
            rad = zeros(K, 1);
            for k = 1:K
                G_k = create_G(t_L, theta_LK(:, k), phi_LK(:, k), M_Lk, N_L, lambda);
                [grad(:, k), bar_g_Lkn] = cal_gradient_f1_t_L(t_L, theta_LK(:, k), phi_LK(:, k), lambda, ...
                                                     z(k), w_L, Sigma_LK(:,:,k), G_k, n);
                [~, tau(k)] = cal_hessian_f1_t_L(t_L, theta_LK(:, k), phi_LK(:, k), lambda, ...
                                                        z(k), w_L, Sigma_LK(:,:,k), bar_g_Lkn, n);
                t_L_k_star(:, k) = t_L(:, n) - real(grad(:, k)) / real(tau(k));
                rad(k) = norm(real(grad(:, k)) / real(tau(k)), 2);
                % G_GK(:, :, k) = G_k;
            end

            G_LA = create_G(t_L, theta_LA, phi_LA, M_LA, N_L, lambda);

            [grad_LA, bar_g_LA] = cal_gradient_f2(t_L, theta_LA, phi_LA, lambda, ...
                                    v_A, w_L, F_LA, Sigma_LA, G_LA, n);
            [~, tau_LA] = cal_hessian_f2(t_L, theta_LA, phi_LA, lambda, ...
                                  v_A, w_L, F_LA, Sigma_LA, bar_g_LA, n);
            % tau_LA = 1e-3;
            t_Lstar = t_L(:, n) - real(grad_LA) / real(tau_LA);
            % temp = t_G(:, n);
            if n==1
                t_L_opt = t_L(:, 2:end);
                t_Ln = find_pos_tL(t_L_opt.', t_L(:, n).', t_Lstar.', lambda, normA, norm(real(grad_LA) / real(tau_LA), 2), t_L_k_star.', rad);
                if isnan(t_Ln(1,1)) || isnan(t_Ln(1,2))
                    t_L(:, n) = t_L(:, n);
                else
                    t_L(:, n) = t_Ln.';
                end
            else
                t_L_opt = t_L(:, 1:n-1);
                if n<N_L
                    t_L_opt = [t_L_opt, t_L(:, n+1:end)];
                end
                t_Ln = find_pos_tL(t_L_opt.', t_L(:, n).', t_Lstar.', lambda, normA, norm(real(grad_LA) / real(tau_LA), 2), t_L_k_star.', rad);
                if isnan(t_Ln(1,1)) || isnan(t_Ln(1,2))
                    t_L(:, n) = t_L(:, n);
                else
                    t_L(:, n) = t_Ln.';
                end
            end
            % f(t_Gn.')
        end
        
    end

end