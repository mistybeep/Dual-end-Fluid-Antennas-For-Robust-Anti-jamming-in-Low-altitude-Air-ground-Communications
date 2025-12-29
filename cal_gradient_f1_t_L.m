function [gradient, bar_g_Lkn] = cal_gradient_f1_t_L(t_L, theta_Lk, phi_Lk, lambda, ...
                                                     z_k, w_L, Sigma_LK, G_LK, n)
% 初始化输出
gradient = zeros(2, 1);

N_L = size(w_L, 1);
M_Lk = size(theta_Lk, 1);
Sigma_LK_bar = ones(1,M_Lk) * Sigma_LK;
Sigma_LK_tilde = Sigma_LK_bar' * Sigma_LK_bar;

% 计算g_dot和g_prime
g_dot = zeros(M_Lk, 1);
g_prime = zeros(M_Lk, 1);
for m = 1:M_Lk
    n_Lk = [cos(theta_Lk(m)) * cos(phi_Lk(m)); sin(theta_Lk(m))];
    phase = 2*pi/lambda * (n_Lk' * t_L(:, n));
    g_dot(m) = cos(theta_Lk(m)) * cos(phi_Lk(m)) * exp(1j * phase);
    g_prime(m) = sin(theta_Lk(m)) * exp(1j * phase);
end

% 计算bar{g}_G,k,n
% bar_g_Gkn = zeros(M_Gk, 1);
sum_w_term = zeros(N_L, 1);
for i = 1:N_L
    sum_w_term(i) = conj(w_L(n)) * w_L(i);
end

temp_sum = 0;
for i = 1:N_L
    temp_sum = temp_sum + sum_w_term(i) * Sigma_LK_tilde * G_LK(:,i);
end
bar_g_Lkn = z_k^2 * temp_sum;

% 计算梯度分量
sum_x = 0;
sum_y = 0;
for m = 1:M_Lk
    sum_x = sum_x + 10*imag(conj(bar_g_Lkn(m)) * g_dot(m));
    sum_y = sum_y + 10*imag(conj(bar_g_Lkn(m)) * g_prime(m));
end

gradient(1) = -4*pi/lambda * sum_x;  % ∂f_1k/∂x_Gn
gradient(2) = -4*pi/lambda * sum_y;  % ∂f_1k/∂y_Gn

end