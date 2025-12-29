function [gradient, bar_g_Gkn] = cal_gradient_f1(t_G, theta_Gk, phi_Gk, lambda, ...
                                 z_k, y_k, W_G, Sigma, G, k, n)
% 初始化输出
gradient = zeros(2, 1);

[N_G, K] = size(W_G);
M_Gk = size(theta_Gk, 1);
w_Gk = W_G(:, k);
Sigma_bar = ones(1,M_Gk) * Sigma;
Sigma_tilde = Sigma_bar' * Sigma_bar;

% 计算g_dot和g_prime
g_dot = zeros(M_Gk, 1);
g_prime = zeros(M_Gk, 1);
for m = 1:M_Gk
    n_Gk = [cos(theta_Gk(m)) * cos(phi_Gk(m)); sin(theta_Gk(m))];
    phase = 2*pi/lambda * (n_Gk' * t_G(:, n));
    g_dot(m) = cos(theta_Gk(m)) * cos(phi_Gk(m)) * exp(1j * phase);
    g_prime(m) = sin(theta_Gk(m)) * exp(1j * phase);
end

% 计算bar{g}_G,k,n
% bar_g_Gkn = zeros(M_Gk, 1);
sum_w_term = zeros(N_G, 1);
for l = 1:K
    w_Gl = W_G(:, l);
    for i = 1:N_G
        sum_w_term(i) = sum_w_term(i) + conj(w_Gl(n)) * w_Gl(i);
    end
end

temp_sum = 0;
for i = 1:N_G
    temp_sum = temp_sum + sum_w_term(i) * Sigma_tilde * G(:,i);
end
bar_g_Gkn = z_k^2 * temp_sum - z_k * sqrt(1 + y_k) * conj(w_Gk(n)) * Sigma_bar';

% 计算梯度分量
sum_x = 0;
sum_y = 0;
for m = 1:M_Gk
    sum_x = sum_x + imag(conj(bar_g_Gkn(m)) * g_dot(m));
    sum_y = sum_y + imag(conj(bar_g_Gkn(m)) * g_prime(m));
end

gradient(1) = -4*pi/lambda * sum_x;  % ∂f_1k/∂x_Gn
gradient(2) = -4*pi/lambda * sum_y;  % ∂f_1k/∂y_Gn

end