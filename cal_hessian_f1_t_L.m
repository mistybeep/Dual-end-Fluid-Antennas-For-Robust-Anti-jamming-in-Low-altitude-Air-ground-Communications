function [hessian, max_eigenvalue] = cal_hessian_f1_t_L(t_L, theta_Lk, phi_Lk, lambda, ...
                                                        z_k, w_L, Sigma_LK, bar_g_Lk, n)

% 初始化Hessian矩阵
hessian = zeros(2, 2);

M_Lk = size(theta_Lk, 1);
Sigma_LK_bar = ones(1,M_Lk) * Sigma_LK;
Sigma_LK_tilde = Sigma_LK_bar' * Sigma_LK_bar;

% 计算常数项
const = 8 * pi^2 / lambda^2;

w_sum = abs(w_L(n))^2;

% 初始化中间变量
g_dot = zeros(M_Lk, 1);      % \dot{g}
g_prime = zeros(M_Lk, 1);    % g'
g_ddot = zeros(M_Lk, 1);     % \ddot{g}
g_dot_prime = zeros(M_Lk, 1); % \dot{g}'
g_dprime = zeros(M_Lk, 1);   % g''

% 计算各个导数项
for m = 1:M_Lk
    n_Lk = [cos(theta_Lk(m)) * cos(phi_Lk(m)); sin(theta_Lk(m))];
    phase = 2*pi/lambda * (n_Lk' * t_L(:, n));
    cos_theta = cos(theta_Lk(m));
    sin_theta = sin(theta_Lk(m));
    cos_phi = cos(phi_Lk(m));
    
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

for m = 1:M_Lk
    % ∂²f/∂x²
    term1_xx = z_k^2 * w_sum * real(conj(Sigma_LK_tilde(m,:)) * conj(g_dot) * g_dot(m));
    term2_xx = real(conj(bar_g_Lk(m)) * g_ddot(m));
    sum_xx = sum_xx + (term1_xx - term2_xx);
    
    % ∂²f/∂x∂y
    term1_xy = z_k^2 * w_sum * real(conj(Sigma_LK_tilde(m,:)) * conj(g_prime) * g_dot(m));
    term2_xy = real(conj(bar_g_Lk(m)) * g_dot_prime(m));
    sum_xy = sum_xy + (term1_xy - term2_xy);
    
    % ∂²f/∂y∂x
    term1_yx = z_k^2 * w_sum * real(conj(Sigma_LK_tilde(m,:)) * conj(g_dot) * g_prime(m));
    term2_yx = real(conj(bar_g_Lk(m)) * g_dot_prime(m));
    sum_yx = sum_yx + (term1_yx - term2_yx);
    
    % ∂²f/∂y²
    term1_yy = z_k^2 * w_sum * real(conj(Sigma_LK_tilde(m,:)) * conj(g_prime) * g_prime(m));
    term2_yy = real(conj(bar_g_Lk(m)) * g_dprime(m));
    sum_yy = sum_yy + (term1_yy - term2_yy);
end

% 填充Hessian矩阵
hessian(1,1) = const * sum_xx;  % ∂²f/∂x²
hessian(1,2) = const * sum_xy;  % ∂²f/∂x∂y
hessian(2,1) = const * sum_yx;  % ∂²f/∂y∂x
hessian(2,2) = const * sum_yy;  % ∂²f/∂y²

% 计算 Hessian 矩阵的特征值
eigenvalues = eig(hessian);

max_eigenvalue = max(eigenvalues);
if max_eigenvalue < 0
    max_eigenvalue = 1e4;
else
    % max_eigenvalue = max(max_eigenvalue, 0) + 1e-4;
    order_of_magnitude = floor(log10(abs(max_eigenvalue) + 1e-10));
    max_eigenvalue = max_eigenvalue + 10^(order_of_magnitude + 1);
end
end