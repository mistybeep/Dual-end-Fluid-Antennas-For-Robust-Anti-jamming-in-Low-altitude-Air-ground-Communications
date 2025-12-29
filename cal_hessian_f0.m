function [hessian, max_eigenvalue] = cal_hessian_f0(t_G, theta_GA, phi_GA, lambda, ...
                                  v_A, W_G, F, Sigma, bar_g_GA, n)

% 初始化Hessian矩阵
hessian = zeros(2, 2);

[~, K] = size(W_G);
M_GA = size(theta_GA, 1);
Sigma_bar = F' * Sigma;
Sigma_tilde = Sigma_bar' * (v_A * v_A') * Sigma_bar;
% 计算常数项
const = 8 * pi^2 / lambda^2;

% 计算权重平方和
w_sum = 0;
for l = 1:K
    w_Gl = W_G(:, l);
    w_sum = w_sum + abs(w_Gl(n))^2;
end

% 初始化中间变量
g_dot = zeros(M_GA, 1);      % \dot{g}
g_prime = zeros(M_GA, 1);    % g'
g_ddot = zeros(M_GA, 1);     % \ddot{g}
g_dot_prime = zeros(M_GA, 1); % \dot{g}'
g_dprime = zeros(M_GA, 1);   % g''

% 计算各个导数项
for m = 1:M_GA
    n_GA = [cos(theta_GA(m)) * cos(phi_GA(m)); sin(theta_GA(m))];
    phase = 2*pi/lambda * (n_GA' * t_G(:, n));
    cos_theta = cos(theta_GA(m));
    sin_theta = sin(theta_GA(m));
    cos_phi = cos(phi_GA(m));
    
    g_dot(m) = cos_theta * cos_phi * exp(1j * phase);
    g_prime(m) = sin_theta * exp(1j * phase);
    g_ddot(m) = cos_theta^2 * cos_phi^2 * exp(1j * phase);
    g_dot_prime(m) = cos_theta * sin_theta * cos_phi * exp(1j * phase);
    g_dprime(m) = sin_theta^2 * exp(1j * phase);
end

% 计算Hessian矩阵的各个元素
sum_xx = 0;
sum_xy = 0;
sum_yx = 0;
sum_yy = 0;

for m = 1:M_GA
    % ∂²f/∂x²
    term1_xx = w_sum * real(conj(Sigma_tilde(m,:)) * conj(g_dot) * g_dot(m));
    term2_xx = real(conj(bar_g_GA(m)) * g_ddot(m));
    sum_xx = sum_xx + (term1_xx - term2_xx);
    
    % ∂²f/∂x∂y
    term1_xy = w_sum * real(conj(Sigma_tilde(m,:)) * conj(g_prime) * g_dot(m));
    term2_xy = real(conj(bar_g_GA(m)) * g_dot_prime(m));
    sum_xy = sum_xy + (term1_xy - term2_xy);
    
    % ∂²f/∂y∂x
    term1_yx = w_sum * real(conj(Sigma_tilde(m,:)) * conj(g_dot) * g_prime(m));
    term2_yx = real(conj(bar_g_GA(m)) * g_dot_prime(m));
    sum_yx = sum_yx + (term1_yx - term2_yx);
    
    % ∂²f/∂y²
    term1_yy = w_sum * real(conj(Sigma_tilde(m,:)) * conj(g_prime) * g_prime(m));
    term2_yy = real(conj(bar_g_GA(m)) * g_dprime(m));
    sum_yy = sum_yy + (term1_yy - term2_yy);
end

% 填充Hessian矩阵
hessian(1,1) = const * sum_xx;  % ∂²f/∂x²
hessian(1,2) = const * sum_xy;  % ∂²f/∂x∂y
hessian(2,1) = const * sum_yx;  % ∂²f/∂y∂x
hessian(2,2) = const * sum_yy;  % ∂²f/∂y²

% 计算 Hessian 矩阵的特征值
eigenvalues = eig(hessian);

% 求最大特征值
max_eigenvalue = max(eigenvalues);
max_eigenvalue = max(max_eigenvalue, 0) + 1e-4;

end