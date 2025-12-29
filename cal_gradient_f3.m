function [gradient, bar_f_A, bar_f_JA] = cal_gradient_f3(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                    theta_LA, phi_LA, Sigma_LA, F_GA, F_LA, ...
                                    theta_JA, phi_JA, Sigma_JA, P_J, ...
                                    v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                                    mu, nu, kappa, n, lambda, err_range)
% 初始化输出
gradient = zeros(2, 1);

N_A = size(v_A, 1);
M_GA = size(theta_GA, 1);
M_LA = size(theta_LA, 1);
M_JA = size(theta_JA, 1);

M_k = [M_GA; M_LA];

Q = size(mu, 1);
R = size(nu, 1);

theta = cell(2,1);
theta{1} = theta_GA; 
theta{2} = theta_LA;

phi = cell(2,1);
phi{1} = phi_GA; 
phi{2} = phi_LA;

% Sigma = cell(2,1);
% Sigma{1} = Sigma_GA; 
% Sigma{2} = Sigma_LA;

Sigma_GA_bar = Sigma_GA * G_GA;
Sigma_GA_tilde = Sigma_GA_bar * (W_G * W_G') * Sigma_GA_bar';

Sigma_LA_bar = Sigma_LA * G_LA;
Sigma_LA_tilde = Sigma_LA_bar * (w_L * w_L') * Sigma_LA_bar';

SigmaJA_bar = Sigma_JA * G_JA;
Sigma_JA_tilde = SigmaJA_bar * SigmaJA_bar';


bar_f_A = cell(2,1);

sum_GA = 0;
sum_LA = 0;
for i = 1:N_A
    sum_GA = sum_GA + kappa * v_A(i) * conj(v_A(n)) * Sigma_GA_tilde * F_GA(:,i);
    sum_LA = sum_LA - v_A(i) * conj(v_A(n)) * Sigma_LA_tilde * F_LA(:,i);
end
bar_f_A{1} = sum_GA;
bar_f_A{2} = sum_LA;

% 计算A部分的梯度
sum_x_A = 0;
sum_y_A = 0;
for k = 1:2
    for m = 1:M_k(k)
        n_A = [cos(theta{k}(m)) * cos(phi{k}(m)); cos(theta{k}(m)) * sin(phi{k}(m))];
        phase = 2*pi/lambda * (n_A' * r_A(:, n));
        f_dot = cos(theta{k}(m)) * cos(phi{k}(m)) * exp(1j * phase);
        f_prime = cos(theta{k}(m)) * sin(phi{k}(m)) * exp(1j * phase);
        
        sum_x_A = sum_x_A + imag(conj(bar_f_A{k}(m)) * f_dot);
        sum_y_A = sum_y_A + imag(conj(bar_f_A{k}(m)) * f_prime);
    end
end

bar_f_JA = zeros(M_JA, Q, R);
sum_x_JA = 0;
sum_y_JA = 0;
for q = 1:Q
    theta_q = theta_JA - err_range/2 + (q-1)*err_range/(Q-1);
    for r = 1:R
        phi_r = phi_JA - err_range/2 + (r-1)*err_range/(R-1);
        F_JA = zeros(M_JA, N_A);
        for m = 1:M_JA
            n_r = [cos(theta_q(m)) * cos(phi_r(m)); cos(theta_q(m)) * sin(phi_r(m))];
            for nr = 1:N_A
                rho_r = n_r' * r_A(:, nr);
                F_JA(m, nr) = exp(1j * 2 * pi * rho_r / lambda);
            end
        end
        sum_JA = 0;
        for i = 1:N_A
            sum_JA = sum_JA + kappa * P_J * v_A(i) * conj(v_A(n)) * ...
                     Sigma_JA_tilde * F_JA(:,i);
        end
        bar_f_JA(:,q,r) = sum_JA;
        for m = 1:M_JA
            n_A = [cos(theta_q(m)) * cos(phi_r(m)); cos(theta_q(m)) * sin(phi_r(m))];
            phase = 2*pi/lambda * (n_A' * r_A(:, n));
            f_dot_JA = cos(theta_q(m)) * cos(phi_r(m)) * exp(1j * phase);
            f_prime_JA = cos(theta_q(m)) * sin(phi_r(m)) * exp(1j * phase);

            sum_x_JA = sum_x_JA + mu(q) * nu(r) * ...
                       imag(conj(bar_f_JA(m,q,r)) * f_dot_JA);
            sum_y_JA = sum_y_JA + mu(q) * nu(r) * ...
                       imag(conj(bar_f_JA(m,q,r)) * f_prime_JA);
        end
    end
end

% 计算JA部分的梯度
% sum_x_JA = 0;
% sum_y_JA = 0;
% for q = 1:Q
%     theta_q = theta_JA - err_range/2 + (q-1)*err_range/(Q-1);
%     for r = 1:R
%         phi_r = phi_JA - err_range/2 + (r-1)*err_range/(R-1);
%         for m = 1:M_JA
%             n_A = [cos(theta_q(m)) * cos(phi_r(m)); cos(theta_q(m)) * sin(phi_r(m))];
%             phase = 2*pi/lambda * (n_A' * r_A(:, n));
%             f_dot_JA = cos(theta_q(m)) * cos(phi_r(m)) * exp(1j * phase);
%             f_prime_JA = cos(theta_q(m)) * sin(phi_r(m)) * exp(1j * phase);
% 
%             sum_x_JA = sum_x_JA + mu(q) * nu(r) * ...
%                        imag(conj(bar_f_JA(m,q,r)) * f_dot_JA);
%             sum_y_JA = sum_y_JA + mu(q) * nu(r) * ...
%                        imag(conj(bar_f_JA(m,q,r)) * f_prime_JA);
%         end
%     end
% end

% 组合总梯度
const = -4*pi/lambda;
gradient(1) = const * (sum_x_A + sum_x_JA);  % ∂f₃/∂x_An
gradient(2) = const * (sum_y_A + sum_y_JA);  % ∂f₃/∂y_An

end