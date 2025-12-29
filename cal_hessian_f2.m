function [hessian, max_eigenvalue] = cal_hessian_f2(t_L, theta_LA, phi_LA, lambda, ...
                                  v_A, w_L, F_LA, Sigma_LA, bar_g_LA, n)

% 初始化Hessian矩阵
hessian = zeros(2, 2);

M_LA = size(theta_LA, 1);
Sigma_LA_bar = F_LA' * Sigma_LA;
Sigma_LA_tilde = Sigma_LA_bar' * (v_A * v_A') * Sigma_LA_bar;
% 计算常数项
const = 8 * pi^2 / lambda^2;

% 计算权重平方和
w_sum = abs(w_L(n))^2;

% 初始化中间变量
g_dot = zeros(M_LA, 1);      % \dot{g}
g_prime = zeros(M_LA, 1);    % g'
g_ddot = zeros(M_LA, 1);     % \ddot{g}
g_dot_prime = zeros(M_LA, 1); % \dot{g}'
g_dprime = zeros(M_LA, 1);   % g''

% 计算各个导数项
for m = 1:M_LA
    n_LA = [cos(theta_LA(m)) * cos(phi_LA(m)); sin(theta_LA(m))];
    phase = 2*pi/lambda * (n_LA' * t_L(:, n));
    cos_theta = cos(theta_LA(m));
    sin_theta = sin(theta_LA(m));
    cos_phi = cos(phi_LA(m));
    
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

for m = 1:M_LA
    % ∂²f/∂x²
    term1_xx = -1 * w_sum * real(conj(Sigma_LA_tilde(m,:)) * conj(g_dot) * g_dot(m));
    term2_xx = real(conj(bar_g_LA(m)) * g_ddot(m));
    sum_xx = sum_xx + (term1_xx - term2_xx);
    
    % ∂²f/∂x∂y
    term1_xy = -1 * w_sum * real(conj(Sigma_LA_tilde(m,:)) * conj(g_prime) * g_dot(m));
    term2_xy = real(conj(bar_g_LA(m)) * g_dot_prime(m));
    sum_xy = sum_xy + (term1_xy - term2_xy);
    
    % ∂²f/∂y∂x
    term1_yx =  -1 * w_sum * real(conj(Sigma_LA_tilde(m,:)) * conj(g_dot) * g_prime(m));
    term2_yx = real(conj(bar_g_LA(m)) * g_dot_prime(m));
    sum_yx = sum_yx + (term1_yx - term2_yx);
    
    % ∂²f/∂y²
    term1_yy =  -1 * w_sum * real(conj(Sigma_LA_tilde(m,:)) * conj(g_prime) * g_prime(m));
    term2_yy = real(conj(bar_g_LA(m)) * g_dprime(m));
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