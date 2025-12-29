function Sigma = generate_Sigma(M, rho_0, alpha, d)
    variance = rho_0 * d^(-alpha);
    Sigma = zeros(M, M);
    diag_elements = sqrt(variance) * (randn(M, 1) + 1j * randn(M, 1)) / sqrt(2); % CSCG 分布
    % diag_elements = (randn(M, 1) + 1j * randn(M, 1)) / sqrt(2);
    Sigma(1, 1) = diag_elements(1); % 第一个元素为 LoS 分量
    for m = 2:M
        Sigma(m, m) = diag_elements(m) / (M - 1); % 非 LoS 分量除以 M-1
    end

end