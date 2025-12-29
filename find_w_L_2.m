function w_L = find_w_L_2(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2)

    tildeR = zeros(K, 1);
    for k = 1 : K
        tildeR(k) = log2(1 + y(k)) - y(k) ...
                    + 2 * z(k) * sqrt(1 + y(k)) * real(H_G(k,:) * W_G(:, k)) ...
                    - z(k)^2 * (norm(H_G(k,:) * W_G, 2)^2 + sigma2) - gamma(k);
    end

    Hbar = H_LA'*(v_A*v_A')*H_LA;
    
    N_L = size(H_L, 2);
    Phi_cell = cell(K,1);
    for k = 1 : K
        Phi_cell{k} = (z(k)^2) * (H_L(k, :)' * H_L(k, :));
    end

    % 收敛参数（可调）
    Imax = 500;
    eps_k = 1e-5;        % epsilon for convergence
    step_a = 1; step_b = 0.5; step_c = 10.0;   % step size params: t_n = c/(a + b*n)

    nu = zeros(K+1,1);    % nu(1)=nu0, nu(2..K+1)=nu1..nuK
    nu(1) = 1;
    n = 0;

    Phi_inv_sqrt = eye(N_L)^(-0.5);
    A = Phi_inv_sqrt * Hbar * Phi_inv_sqrt';

    [U_A, D_A] = eig(A);
    [~, idx_max] = max(real(diag(D_A)));
    u = U_A(:, idx_max);  % 主特征向量（列）
    u = u / norm(u);      % 单位化
    w_L = sqrt(P) * Phi_inv_sqrt * u;
    
    overlineP = P + sum(tildeR);
    flag = 1;
    for k = 1 : K
        if real(trace(Phi_cell{k} * (w_L * w_L'))) > tildeR(k)
            flag = 0;
        end
    end
    if flag == 0
        while true
            % --- 构造 Phi = nu0 * I + sum_{i=1..K} nu_i * Phi_i
            beta = nu / (nu(1) * P + sum(nu(2:K+1).* tildeR));
            Phi = beta(1) * eye(N_L);
            for k = 1:K
                Phi = Phi + beta(1+k) * Phi_cell{k};
            end
            Phi_inv_sqrt = Phi^(-0.5);
            A = Phi_inv_sqrt * Hbar * Phi_inv_sqrt;
            % 确保 A 是 Hermitian（数值误差）
            A = (A + A')/2;

            % --- 主特征向量 u (largest eigenvalue)
            [U_A, D_A] = eig(A);
            [~, idx_max] = max(real(diag(D_A)));
            u = U_A(:, idx_max);  % 主特征向量（列）
            u = u / norm(u);      % 单位化
            w_L = sqrt(overlineP) * Phi_inv_sqrt * u;

            % --- 还原回 W_L = Phi^{-1/2} * W_tilde * Phi^{-1/2}
            WL = w_L * w_L';
            % --- 计算次梯度 g: g0 = Tr(I*WL) - P_L_max; gi = Tr(Phi_i * WL) - tildeR_i
            g = zeros(K+1,1);
            g(1) = real(trace(eye(N_L) * WL)) - P;
            for k = 1:K
                g(1+k) = real(trace(Phi_cell{k} * WL)) - tildeR(k);
            end

            % --- 终止判定： nu_i * g_i <= eps (按元素)
            cond_vec = nu .* g;
            if all(abs(cond_vec) <= eps_k) || n >= Imax
                break;
            end

            % --- 更新步长 t_n 并更新 nu = [nu + t_n * g]^+
            tn = step_c / (step_a + step_b * n);
            nu = max(nu + tn * g, 1e-6);

            % --- 迭代计数
            n = n + 1;
        end
    end

end