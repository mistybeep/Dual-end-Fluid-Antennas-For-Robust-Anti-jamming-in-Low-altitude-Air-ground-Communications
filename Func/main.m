M = 8;
K = 4;
data = load_dataset('dataset.mat');
alpha = 2.8;
rho_0_dB = -40;
rho_0 = 10^(rho_0_dB / 10); % 转换为线性单位
Sigma_GA_Test = zeros(100,M,M);
Sigma_LA_Test = zeros(100,M,M);
Sigma_JA_Test = zeros(100,M,M);
Sigma_GK_Test = zeros(100,M,M,K);
Sigma_LK_Test = zeros(100,M,M,K);
for i = 1:100
    Sigma_GA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_GA);
    Sigma_LA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_LA);
    Sigma_JA_Test(i,:,:) = generate_Sigma(M, rho_0, alpha, data.d_JA);
    for k = 1:K
        Sigma_GK_Test(i,:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_G(k));
        Sigma_LK_Test(i,:,:,k) = generate_Sigma(M, rho_0, alpha, data.d_single_L(k));
    end
end
save("Sigma.mat", 'Sigma_GA_Test', 'Sigma_LA_Test', 'Sigma_JA_Test', 'Sigma_GK_Test','Sigma_LK_Test')