function [gradient, bar_g_LA] = cal_gradient_f2(t_L, theta_LA, phi_LA, lambda, ...
                                    v_A, w_L, F_LA, Sigma_LA, G_LA, n)
% 初始化输出
gradient = zeros(2, 1);

N_L = size(w_L, 1);
M_LA = size(theta_LA, 1);
Sigma_LA_bar = F_LA' * Sigma_LA;
Sigma_LA_tilde = Sigma_LA_bar' * (v_A * v_A') * Sigma_LA_bar;

g_dot = zeros(M_LA, 1);
g_prime = zeros(M_LA, 1);
for m = 1:M_LA
    n_LA = [cos(theta_LA(m)) * cos(phi_LA(m)); sin(theta_LA(m))];
    phase = 2*pi/lambda * (n_LA' * t_L(:, n));
    g_dot(m) = cos(theta_LA(m)) * cos(phi_LA(m)) * exp(1j * phase);
    g_prime(m) = sin(theta_LA(m)) * exp(1j * phase);
end

% 计算bar{g}_G,k,n
sum_w_term = zeros(N_L, 1);
for i = 1:N_L
    sum_w_term(i) = conj(w_L(n)) * w_L(i);
end


temp_sum = 0;
for i = 1:N_L
    temp_sum = temp_sum - sum_w_term(i) * Sigma_LA_tilde * G_LA(:,i);
end
bar_g_LA = temp_sum;

% 计算梯度分量
sum_x = 0;
sum_y = 0;
for m = 1:M_LA
    sum_x = sum_x + imag(conj(bar_g_LA(m)) * g_dot(m));
    sum_y = sum_y + imag(conj(bar_g_LA(m)) * g_prime(m));
end

gradient(1) = -4*pi/lambda * sum_x;
gradient(2) = -4*pi/lambda * sum_y;

end