function [hessian, max_eigenvalue] = cal_hessian_f3(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                    theta_LA, phi_LA, Sigma_LA, ...
                                    theta_JA, phi_JA, Sigma_JA, P_J, ...
                                    v_A, W_G, w_L, G_GA, G_LA, G_JA, bar_f_A, bar_f_JA, ...
                                    mu, nu, kappa, n, lambda, err_range)

% 初始化Hessian矩阵
hessian = zeros(2, 2);

% N_A = size(v_A, 1);
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

Sigma_GA_bar = Sigma_GA * G_GA;
Sigma_GA_tilde = Sigma_GA_bar * (W_G * W_G') * Sigma_GA_bar';

Sigma_LA_bar = Sigma_LA * G_LA;
Sigma_LA_tilde = Sigma_LA_bar * (w_L * w_L') * Sigma_LA_bar';

SigmaJA_bar = Sigma_JA * G_JA;
Sigma_JA_tilde = SigmaJA_bar * SigmaJA_bar';

Sigma = cell(2,1);
Sigma{1} = Sigma_GA_tilde; 
Sigma{2} = Sigma_LA_tilde;

% 计算常数项
const = 8 * pi^2 / lambda^2;
v_An_sq = abs(v_A(n))^2;

bar_kappa = [kappa; -1];

% 计算A部分的中间变量
f_dot_A = cell(2, 1);
f_prime_A = cell(2, 1);
f_ddot_A = cell(2, 1);
f_dot_prime_A = cell(2, 1);
f_dprime_A = cell(2, 1);

for k = 1:2
    for m = 1:M_k(k)
        n_A = [cos(theta{k}(m)) * cos(phi{k}(m)); cos(theta{k}(m)) * sin(phi{k}(m))];
        phase = 2*pi/lambda * (n_A' * r_A(:, n));
        cos_theta = cos(theta{k}(m));
        sin_theta = sin(theta{k}(m));
        cos_phi = cos(phi{k}(m));
        sin_phi = sin(phi{k}(m));
        
        f_dot_A{k}(m,1) = cos_theta * cos_phi * exp(1j * phase);
        f_prime_A{k}(m,1) = cos_theta * sin_phi * exp(1j * phase);
        f_ddot_A{k}(m,1) = cos_theta^2 * cos_phi^2 * exp(1j * phase);
        f_dot_prime_A{k}(m,1) = cos_theta^2 * cos_phi * sin_phi * exp(1j * phase);
        f_dprime_A{k}(m,1) = cos_theta^2 * sin_phi^2 * exp(1j * phase);
    end
end

f_dot_JA = zeros(M_JA, Q, R);
f_prime_JA = zeros(M_JA, Q, R);
f_ddot_JA = zeros(M_JA, Q, R);
f_dot_prime_JA = zeros(M_JA, Q, R);
f_dprime_JA = zeros(M_JA, Q, R);

for q = 1:Q
    theta_q = theta_JA - err_range/2 + (q-1)*err_range/(Q-1);
    for r = 1:R
        phi_r = phi_JA - err_range/2 + (r-1)*err_range/(R-1);
        for m = 1:M_JA
            n_A = [cos(theta_q(m)) * cos(phi_r(m)); cos(theta_q(m)) * sin(phi_r(m))];
            phase = 2*pi/lambda * (n_A' * r_A(:, n));
            cos_theta = cos(theta_q(m));
            cos_phi = cos(phi_r(m));
            sin_phi = sin(phi_r(m));
            
            f_dot_JA(m,q,r) = cos_theta * cos_phi * exp(1j * phase);
            f_prime_JA(m,q,r) = cos_theta * sin_phi * exp(1j * phase);
            f_ddot_JA(m,q,r) = cos_theta^2 * cos_phi^2 * exp(1j * phase);
            f_dot_prime_JA(m,q,r) = cos_theta^2 * cos_phi * sin_phi * exp(1j * phase);
            f_dprime_JA(m,q,r) = cos_theta^2 * sin_phi^2 * exp(1j * phase);
        end
    end
end

sum_xx = 0;
sum_xy = 0;
sum_yx = 0;
sum_yy = 0;

% A部分贡献
for k = 1:2
    for m = 1:M_k(k)
        % ∂²f/∂x²
        term1_xx = bar_kappa(k) * v_An_sq * real(conj(Sigma{k}(m,:)) * ...
                   conj(f_dot_A{k}) * f_dot_A{k}(m));
        term2_xx = real(conj(bar_f_A{k}(m)) * f_ddot_A{k}(m));
        sum_xx = sum_xx + (term1_xx - term2_xx);
        
        % ∂²f/∂x∂y
        term1_xy = bar_kappa(k) * v_An_sq * real(conj(Sigma{k}(m,:)) * ...
                   conj(f_prime_A{k}) * f_dot_A{k}(m));
        term2_xy = real(conj(bar_f_A{k}(m)) * f_dot_prime_A{k}(m));
        sum_xy = sum_xy + (term1_xy - term2_xy);
        
        % ∂²f/∂y∂x
        term1_yx = bar_kappa(k) * v_An_sq * real(conj(Sigma{k}(m,:)) * ...
                   conj(f_dot_A{k}) * f_prime_A{k}(m));
        term2_yx = real(conj(bar_f_A{k}(m)) * f_dot_prime_A{k}(m));
        sum_yx = sum_yx + (term1_yx - term2_yx);
        
        % ∂²f/∂y²
        term1_yy = bar_kappa(k) * v_An_sq * real(conj(Sigma{k}(m,:)) * ...
                   conj(f_prime_A{k}) * f_prime_A{k}(m));
        term2_yy = real(conj(bar_f_A{k}(m)) * f_dprime_A{k}(m));
        sum_yy = sum_yy + (term1_yy - term2_yy);
    end
end

% JA部分贡献
for q = 1:Q
    for r = 1:R
        for m = 1:M_JA
                % ∂²f/∂x²
            term1_xx = mu(q) * nu(r) * kappa * P_J * v_An_sq * real(conj(Sigma_JA_tilde(m,:)) * ...
                conj(f_dot_JA(:,q,r)) * f_dot_JA(m,q,r));
            term2_xx = mu(q) * nu(r) *real(conj(bar_f_JA(m,q,r)) * f_ddot_JA(m,q,r));
            sum_xx = sum_xx + (term1_xx - term2_xx);

            % ∂²f/∂x∂y
            term1_xy = mu(q) * nu(r) * kappa * P_J * v_An_sq * real(conj(Sigma_JA_tilde(m,:)) * ...
                conj(f_prime_JA(:,q,r)) * f_dot_JA(m,q,r));
            term2_xy = mu(q) * nu(r) *real(conj(bar_f_JA(m,q,r)) * f_dot_prime_JA(m,q,r));
            sum_xy = sum_xy + (term1_xy - term2_xy);

            % ∂²f/∂y∂x
            term1_yx = mu(q) * nu(r) * kappa * P_J * v_An_sq * real(conj(Sigma_JA_tilde(m,:)) * ...
                conj(f_dot_JA(:,q,r)) * f_prime_JA(m,q,r));
            term2_yx = mu(q) * nu(r) *real(conj(bar_f_JA(m,q,r)) * f_dot_prime_JA(m,q,r));
            sum_yx = sum_yx + (term1_yx - term2_yx);

            % ∂²f/∂y²
            term1_yy = mu(q) * nu(r) * kappa * P_J * v_An_sq * real(conj(Sigma_JA_tilde(m,:)) * ...
                conj(f_prime_JA(:,q,r)) * f_prime_JA(m,q,r));
            term2_yy = mu(q) * nu(r) *real(conj(bar_f_JA(m,q,r)) * f_dprime_JA(m,q,r));
            sum_yy = sum_yy + (term1_yy - term2_yy);
        end
    end
end

% 填充Hessian矩阵
hessian(1,1) = const * sum_xx;  % ∂²f/∂x²
hessian(1,2) = const * sum_xy;  % ∂²f/∂x∂y
hessian(2,1) = const * sum_yx;  % ∂²f/∂y∂x
hessian(2,2) = const * sum_yy;  % ∂²f/∂y²

% 计算 Hessian 矩阵的特征值
eigenvalues = real(eig(hessian));

% 求最大特征值
max_eigenvalue = max(eigenvalues);
if max_eigenvalue < 0
    max_eigenvalue = 1e4;
else
    % max_eigenvalue = max(max_eigenvalue, 0) + 1e-4;
    order_of_magnitude = floor(log10(abs(max_eigenvalue) + 1e-10));
    max_eigenvalue = max_eigenvalue + 10^order_of_magnitude;
end
% max_eigenvalue = max(max_eigenvalue, 0) + 1e-5;

end