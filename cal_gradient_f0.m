function [gradient, bar_g_GA] = cal_gradient_f0(t_G, theta_GA, phi_GA, lambda, ...
                                    v_A, W_G, F, Sigma, G, n)
% 初始化输出
gradient = zeros(2, 1);

[N_G, K] = size(W_G);
M_GA = size(theta_GA, 1);
Sigma_bar = F' * Sigma;
Sigma_tilde = Sigma_bar' * (v_A * v_A') * Sigma_bar;

% 计算g_dot和g_prime
g_dot = zeros(M_GA, 1);
g_prime = zeros(M_GA, 1);
for m = 1:M_GA
    n_GA = [cos(theta_GA(m)) * cos(phi_GA(m)); sin(theta_GA(m))];
    phase = 2*pi/lambda * (n_GA' * t_G(:, n));
    g_dot(m) = cos(theta_GA(m)) * cos(phi_GA(m)) * exp(1j * phase);
    g_prime(m) = sin(theta_GA(m)) * exp(1j * phase);
end

% 计算bar{g}_G,k,n
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
bar_g_GA = temp_sum;

% 计算梯度分量
sum_x = 0;
sum_y = 0;
for m = 1:M_GA
    sum_x = sum_x + imag(conj(bar_g_GA(m)) * g_dot(m));
    sum_y = sum_y + imag(conj(bar_g_GA(m)) * g_prime(m));
end

gradient(1) = -4*pi/lambda * sum_x;
gradient(2) = -4*pi/lambda * sum_y;

end