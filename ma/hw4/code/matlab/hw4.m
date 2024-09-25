%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses the LiveLink feature therefore you need to have a
% running "COMSOL Multiphysics 6.0 with MATLAB" instance. Setup guide:
% https://doc.comsol.com/6.0/doc/com.comsol.help.llmatlab/LiveLinkForMATLABUsersGuide.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

addpath('../comsol');

%% Load parameters and impedance data from COMSOL
model = mphopen('soundboard.mph');
Z_table = mphtable(model, 'tbl1');
Z_data = Z_table.data;

Lx = mphevaluate(model, 'Lx');
Ly = mphevaluate(model, 'Ly');
f_step = mphevaluate(model, 'fstep');

%% b) Soundboard design

str_freqs_start = [349.23; 440; 523.25; 659.25; 783.99];
str_freqs = f_step*round(str_freqs_start(:)/f_step);
freq_colors = [
    1 0 0; ... % red
    0 1 0; ... % green
    0 0 1; ... % blue
    0 1 1; ... % cyan
    1 0 1; ... % magenta
    1 1 0; ... % yellow
    0 0 0 ...  % black
];

% Build whole soundboard impedance exploiting simmetry
min_points = 15;
soundboard = [];
for k = 1:length(str_freqs)
    curr_freq = str_freqs(k);
    curr_color = freq_colors(k, :);
    freq_filter = Z_data(:, 3)==curr_freq;
    Z_freq_data = Z_data(freq_filter, :);
    [Z_mink, Z_mink_idx] = mink(Z_freq_data(:, 4), min_points);
    Z = Z_freq_data(Z_mink_idx, :);
    Z_x = Z(:, 1);
    Z_y = Z(:, 2);
    Z_freq = Z(:, 3);
    Z_color = repmat(curr_color, min_points, 1);
    Z_values = Z(:, 4);
    Y_values = 1./Z(:, 4);
    for i = 1:min_points
        curr_x = Z_x(i);
        curr_y = Z_y(i);
        curr_Z_value = Z_values(i);
        curr_Y_value = Y_values(i);
        new_x = [Lx - curr_x; curr_x; Lx - curr_x];
        new_y = [curr_y; Ly - curr_y; Ly - curr_y];
        new_freqs = repmat(curr_freq, 3, 1);
        new_colors = repmat(curr_color, 3, 1);
        new_Z_values = repmat(curr_Z_value, 3, 1);
        new_Y_values = repmat(curr_Y_value, 3, 1);
        Z_x = [Z_x; new_x];
        Z_y = [Z_y; new_y];
        Z_freq = [Z_freq; new_freqs];
        Z_color = [Z_color; new_colors];
        Z_values = [Z_values; new_Z_values];
        Y_values = [Y_values; new_Y_values];
    end
    soundboard = [
        soundboard; ...
        Z_x, Z_y, Z_freq, Z_values, Y_values, Z_color
    ];
end

% Plot soundboard impedance
figure;
hold on;
x = soundboard(:, 1);
y = soundboard(:, 2);
values = db(abs(soundboard(:, 4)));
colors = [soundboard(:, 6), soundboard(:, 7), soundboard(:, 8)];
scatter3(x, y, values, 30, colors, 'filled', ...
    'XJitter', 'rand', 'XJitterWidth', 0.01, ...
    'YJitter', 'rand', 'YJitterWidth', 0.01);
xlim([0 Lx]);
ylim([0 Ly]);
view(2);
for k = 1:length(str_freqs)
  H(k) = scatter3(NaN, NaN, NaN, 1, freq_colors(k, :), 'filled');
end
legend(H, {num2str(str_freqs, 'f_0 = %-d Hz')});

% Select bridge points and extract them from soundboard
bridge_points = [
    0.74 0.46 350.0;
    0.82 0.54 440.0;
    0.90 0.62 525.0;
    0.98 0.70 660.0;
    1.06 0.78 785.0
];

bridge = zeros(length(str_freqs), size(soundboard, 2));
d = eps();
for k = 1:length(str_freqs)
    point_idx = all(abs(soundboard(:, 1:3)-[bridge_points(k, 1) ...
        bridge_points(k, 2) bridge_points(k, 3)])<d, 2);
    bridge(k, :) = soundboard(point_idx, :);
end

%% c) Compute the eigenfrequencies of the two strings in the pairs

% Compute eigenvalues
T = 900; % tension [N]
rho = 10.8e-3; % linear density [Kg/m]
Z_0 = sqrt(rho*T); % string characteristic impedance
Y_b = bridge(:, 5);
csi = (1i*Z_0.*Y_b)./pi;
epsilon = linspace(-0.2, 0.2, 1000);
A = zeros(length(str_freqs), length(epsilon), 3);
for i = 1:length(epsilon)
    curr_eps = epsilon(i);
    mu = sqrt(csi.^2 + curr_eps^2);
    a1 = csi + curr_eps + mu;
    a2 = csi + curr_eps - mu;
    A(:, i, 1:3) = [a1, a2, mu];
end

% Compute eigenfrequencies
eig_freqs_1 = str_freqs.*(1 + real(A(:, :, 1)));
eig_freqs_2 = str_freqs.*(1 + real(A(:, :, 2)));

% Find best epsilon value
% TODO - The condition of minimum decay times of a fixed value is missing
% [min_freq_dist, best_eps_idx] = min(abs(eig_freqs_1 - eig_freqs_2), [], 2);
% best_eps = epsilon(1, best_eps_idx);
eps = [0.0046, 0.0046, 0.0046, 0.0017, 0.0046];
fprintf('-----------------------------------------------------------\n');
fprintf('Epsilon\n');
fprintf('-----------------------------------------------------------\n');
fprintf('%.4f\n', eps);
fprintf('-----------------------------------------------------------\n');

% Plot eigenvalues
for i = 1:length(str_freqs)
    curr_a_1 = A(i, :, 1);
    curr_a_2 = A(i, :, 2);
    % Plot eigenvalues real part
    figure('name', sprintf('f = %d', str_freqs(i)));
    subplot(1, 2, 1);
    plot(epsilon, real(curr_a_1));
    hold on;
    plot(epsilon, real(curr_a_2));
    xlabel('\epsilon');
    ylabel('\omega_0 x Re\{a\}');
    % Plot eigenvalues imaginary part
    subplot(1, 2, 2);
    plot(epsilon, imag(curr_a_1));
    hold on;
    plot(epsilon, imag(curr_a_2));
    xlabel('\epsilon');
    ylabel('Im\{a\}');
    set(gca, 'FontSize', 10);
end

% Compute optimum eigenfrequencies
fprintf('Eigenfrequencies\n');
fprintf('-----------------------------------------------------------\n');
eig_freqs = zeros(length(str_freqs), 2);
eig_values = zeros(length(str_freqs), 2);
for i = 1:length(str_freqs)
    curr_freq = str_freqs_start(i);
    curr_eps = eps(i);
    mu = sqrt(csi.^2 + curr_eps^2);
    eig_values(i, 1) = csi(i) + curr_eps + mu(i);
    eig_values(i, 2) = csi(i) + curr_eps - mu(i);
    eig_freqs(i, 1) = str_freqs_start(i)*(1 + real(eig_values(i,1)));
    eig_freqs(i, 2) = str_freqs_start(i)*(1 + real(eig_values(i,2)));
    fprintf('f_0 = %.2f Hz, f_1 = %.2f Hz -> ef_0 = %.2f Hz, ef_1 = %.2f Hz\n', ...
        curr_freq, curr_freq*(1+2*curr_eps), eig_freqs(i, 1), eig_freqs(i, 2));
end
fprintf('-----------------------------------------------------------\n');

%% d) Compute the decay times for the two eigenfrequencies for all the 
% pairs

% Compute velocity profile
F_0 = 1;
N = 2^8;
t = linspace(0, 5, N);
V_b = zeros(length(str_freqs), N);
for i = 1:length(str_freqs)
    str_w = 2*pi*str_freqs(i);
    str_eps = eps(i);
    csi_str = csi(i);
    mu_str = sqrt(csi_str^2 + str_eps^2);
    const = 2*pi*F_0*csi_str/(mu_str*Z_0);
    exp_1 = exp(1i*(str_eps+csi_str)*str_w.*t);
    sin_sum = mu_str.*cos(mu_str*str_w.*t)+1i*csi_str.*sin(mu_str*str_w.*t);
    exp_2 = exp(1i*str_w.*t);
    V_b(i, :) = const.*exp_1.*sin_sum.*exp_2;
    figure();
    plot(t, db(abs(V_b(i, :))));
    xlabel('Time [s]');
    ylabel('V_B [dB]');
    title(sprintf('f_0 = %g [Hz]', str_freqs_start(i)));
end

% Compute T60
fprintf('Decay Time\n');
fprintf('-----------------------------------------------------------\n');
t1 = 0;
t2 = 0.01;
t3 = 2;
t4 = 3;
Vb_ref = -60;
T60 = zeros(length(str_freqs), 2);
for i = 1:length(str_freqs)
    V_b1(i, :) = interp1(t, db(abs(V_b(i, :))), t1);
    V_b2(i, :) = interp1(t, db(abs(V_b(i, :))), t2);
    V_b3(i, :) = interp1(t, db(abs(V_b(i, :))), t3);
    V_b4(i, :) = interp1(t, db(abs(V_b(i, :))), t4); 
    m1(i) = abs(V_b1(i) - V_b2(i))/(t2 - t1);
    m2(i) = abs(V_b3(i) - V_b4(i))/(t4 - t3);
    T60(i, 1) = abs(Vb_ref)/m1(i);
    T60(i, 2) = abs(Vb_ref)/m2(i);
    fprintf('A^T_60 = %.3f s, B^T_60 = %.3f s\n', T60(i, 1), T60(i, 2));
end
fprintf('-----------------------------------------------------------\n');
