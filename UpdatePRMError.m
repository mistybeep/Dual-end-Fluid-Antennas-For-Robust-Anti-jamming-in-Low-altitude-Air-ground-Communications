function DeltaPRM = UpdatePRMError(data, v_A, r_A, Sigma_JA, G_JA, mu, nu, lambda, err_range, DeltaPRM, EpsilonPRM)
    
    N_A = size(v_A, 1);
    L = size(data.theta_JA_r, 1);

    Q = size(mu, 1);

    maxIter = 20;         % 最大迭代次数

    % DeltaPRM = zeros(L, L);
    rho = 0.5;            % 初始惩罚系数
    rhoFactor = 1.6;      % 惩罚系数增长因子
    Jval = zeros(maxIter,1);

    Proj = @(X, EpsilonPRM) X * min(1, EpsilonPRM/norm(X, 'fro'));

    FJA = cell(Q^2,1);
    indexSample = 1;
    for q = 1:Q
        theta_q = data.theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
        for r = 1:Q
            phi_r = data.phi_JA_r - err_range/2 + (r-1)*err_range/(Q-1);
            FJA{indexSample} = create_F(r_A, theta_q, phi_r, L, N_A, lambda);
            indexSample = indexSample + 1;
        end
    end

    %% 迭代更新
    for t = 1:maxIter
        % 计算梯度
        grad = zeros(L, L);
        indexSample = 1;
        for q = 1 : Q
            for r = 1 : Q
                Fp = FJA{indexSample};
                grad = grad + mu(q)*nu(r)*Fp * (v_A*v_A') * Fp' * (Sigma_JA + DeltaPRM) * (G_JA*G_JA');
                indexSample = indexSample + 1;
            end
        end
    
        % 更新 DeltaSigma
        DeltaPRM = Proj(diag(diag(DeltaPRM + grad/rho)), EpsilonPRM);
    
        % 更新惩罚参数
        rho = rho * rhoFactor;
    
        % 计算 jamming power
        J = 0;
        indexSample = 1;
        for q = 1 : Q
            for r = 1 : Q
                Fp = FJA{indexSample};
                J = J + mu(q)*nu(r)*norm(v_A' * Fp' * (Sigma_JA + DeltaPRM) * G_JA)^2;
                indexSample = indexSample + 1;
            end
        end
        Jval(t) = J;
    end
    %% 绘制收敛曲线
    % figure;
    % plot(1:maxIter, Jval,'-o','LineWidth',1.5);
    % xlabel('Iteration');
    % ylabel('Jamming Power');
    % title('Convergence of \Delta\Sigma iterations');
    % grid on;

end