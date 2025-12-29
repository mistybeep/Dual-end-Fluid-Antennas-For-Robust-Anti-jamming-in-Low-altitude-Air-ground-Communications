% 主脚本
rng(42)
clear all;
close all;
addpath("")
num_simulations = 100; % 仿真次数
Error_values = 0:4:8; % P_dBm从10到30，间隔2dB
num_Error_values = length(Error_values);

RatesPerMA = zeros(num_simulations, num_Error_values, 11);
RatesMA = zeros(num_simulations, num_Error_values, 11);
RatesRecMA = zeros(num_simulations, num_Error_values, 11);
RatesTransMA = zeros(num_simulations, num_Error_values, 11);
RatesFPA = zeros(num_simulations, num_Error_values, 11);

lambda = 0.1;
norm_A = 3;
M = 8;
N_G = 16;
N_L = 16;
N_J = 16;
N_A = 4;
K = 4;
P_dBm = 10;
P_J_dBm = 30;
sigma2_dBm = - 80;
gamma = 1 * ones(K, 1);
err_range = 4 / 180 * pi;
alpha = 2.8;
rho_0_dB = -40;
Q = 6;
R = 6;
data = load_dataset('dataset.mat');
EpsKappa = 0.05;

w_L = randn(N_L, 1) + 1j * randn(N_L, 1);
w_L = w_L * sqrt(10^(P_dBm / 10) * 1e-3) / norm(w_L, 2);
v_A = randn(N_A, 1) + 1j*randn(N_A, 1);
v_A = v_A / norm(v_A, 2);

if isempty(gcp('nocreate'))
    numWorkers = 20;
    parpool('local', numWorkers);
    % addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});
else
    fprintf('并行池已存在，直接使用\n');
end

% 运行仿真
tic
for idx = 1:num_Error_values
    err_range = Error_values(idx) / 180 * pi;
    fprintf('当前error = %d (%d/%d)\n', err_range, idx, num_Error_values);
    
    parfor sim = 1:num_simulations
        fprintf('  运行仿真 %d/%d...\n', sim, num_simulations);

        ResultPerMA = PerfectDualEndMA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa);
        RatesPerMA(sim, idx, :) = ResultPerMA;

        ResultMA = DualEndMA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa);
        RatesMA(sim, idx, :) = ResultMA;

        ResultRecMA = RecMA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa);
        RatesRecMA(sim, idx, :) = ResultRecMA;

        ResultTransMA = TransMA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa);
        RatesTransMA(sim, idx, :) = ResultTransMA;

        ResultFPA = FPA(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R, w_L, v_A, EpsKappa);
        RatesFPA(sim, idx, :) = ResultFPA;
    end
end
toc