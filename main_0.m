tic

user_X = 500; % Distance between base station and user
ExNumber = 1000; % Number of iterations

M = 8; % Number of base station antennas

Ps_max_dB = [-5:5:20]; % Maximum transmit power in dB
Ps_max = 10.^(Ps_max_dB/10); % Maximum transmit power in linear scale

BW = 100e6; % Bandwidth in Hz
sigma2 = 10^(-10); % User-induced thermal noise power
sigmar2 = 10^(-10); % Active RIS-induced thermal noise power

f_c = 28e9; % Carrier frequency

K = 6; % Number of users
N = 256; % Number of RIS elements
eta_k = ones(K,1); % User weights

large_fading_AI = 3.5;
large_fading_DI = 3.5;

% Initialize variables
Rsum = zeros(length(Ps_max_dB), ExNumber);
Rsum_noRIS = zeros(length(Ps_max_dB), ExNumber);
Rsum_passive = zeros(length(Ps_max_dB), ExNumber);
Rsum_random = zeros(length(Ps_max_dB), ExNumber);

for a = 1:length(Ps_max_dB)
    fprintf('Iteration %d\n', a);
    parfor b = 1:ExNumber
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Position_generate_2(K, user_X);
        [h_k, f_k, G] = Channel_generate2(K, N, M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c);

        Theta = diag(exp(1j * 2 * pi * rand(N, 1)));
        W = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max(a) / K / M);

        [W, Rsum_noRIS(a, b)] = NoRIS_precoding(M, K, N, Ps_max(a), sigma2, eta_k, W, h_k, f_k, G);
        [W, ~, Rsum_random(a, b)] = random_RIS_precoding(M, K, N, Ps_max(a), sigma2, eta_k, Theta, W, h_k, f_k, G);
        [W, Theta, Rsum_passive(a, b)] = passive_RIS_precoding(M, K, N, Ps_max(a), sigma2, eta_k, Theta, W, h_k, f_k, G);
        Theta = 500 * Theta;
        [W, Theta, Rsum(a, b)] = active_RIS_precoding(M, K, N, Ps_max(a)/2, Ps_max(a)/2, sigma2, sigmar2, eta_k, Theta, W, h_k, f_k, G);
    end
end

Rsum_mean = mean(Rsum, 2);
Rsum_noRIS_mean = mean(Rsum_noRIS, 2);
Rsum_passive_mean = mean(Rsum_passive, 2);
Rsum_random_mean = mean(Rsum_random, 2);

figure;
hold on;
box on;
grid on;
plot(Ps_max_dB, Rsum_mean, '-r^', 'LineWidth', 1.5);
plot(Ps_max_dB, Rsum_passive_mean, '-bo', 'LineWidth', 1.5);
plot(Ps_max_dB, Rsum_random_mean, '-m^', 'LineWidth', 1.5);
plot(Ps_max_dB, Rsum_noRIS_mean, '--k', 'LineWidth', 1.5);

xlabel('Total transmit power $P^{\rm max}$ (dBm)', 'Interpreter', 'latex');
ylabel('Sum-rate (bps/Hz)', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times', 'FontSize', 12);

legend('Active RIS (Algorithm 1)', 'Passive RIS', 'Random phase shift', 'Without RIS', 'Interpreter', 'latex', 'FontSize', 12);
save('main_2_2.mat', 'Rsum_mean', 'Rsum_passive_mean', 'Rsum_random_mean', 'Rsum_noRIS_mean');

toc