clc;
close all;
clearvars;

addpath('./functions');

%% Parameters

% Air density
abs_p = 101.325e3; % at sea level 1 atm = 101.325 kPa
abs_temp = 293.15; % 20 °C = 293.15K
dry_air_mol_mass = 4.81e-26;
k_boltz = physconst('Boltzmann');
rho = (abs_p*dry_air_mol_mass)/(k_boltz*abs_temp);

% Speed of sound
air_bulk_modulus = 1.4e5; % 1.4*10^5 Pa
c = sqrt(air_bulk_modulus/rho);

% Sampling
N = 2^16;
f_min = 0;
f_max = 2000;
f_axis = f_min:f_max/N:f_max-f_max/N;
w_min = 2*pi*f_min;
w_max = 2*pi*f_max;
w = 2*pi*f_axis;
k = w./c;

% Exponential horn parameters
l_exp_horn = 0.35;
m = 4.2;
a_0 = 0.008;
a = a_0*exp(m*l_exp_horn);
S_t = pi*a_0^2;
S_m = pi*a^2;
Z_0 = rho*c/S_t;
% pressure waves does not propagate if k < m! we have an evanescent wave
cut_off_freq_idx = find(k < m, 1, 'last');
cut_off_freq = cut_off_freq_idx * f_max/N;

% Compound horn parameters
l_pipe = 0.85;

%% Impedances

% Conical approximated horn (without external load)
n_sections = linspace(1, 20, 20);
delta_axis = l_exp_horn./n_sections;
Z_1 = zeros(length(n_sections), N);
for i = 1:length(n_sections)
    Z_1(i, :) = approximate_exponential_horn(l_exp_horn, a_0, m, ...
        delta_axis(i), k, 0, rho, c);
end

% Analytical exponential horn
b = sqrt(k.^2 - m^2);
theta = atan(m./b);
Z_2 = Z_0.*1i.*sin(b.*l_exp_horn)./cos(b.*l_exp_horn - theta);

%% a) Evaluate the error e1 as a function of the length δ of the conical 
% sections and plot it in Matlab.

e_1 = zeros(1, length(n_sections));
for i = 1:length(e_1)
    e_1(1, i) = mean(abs(Z_1(i, :)-Z_2).^2, 'omitnan');
end
figure();
stem(flip(delta_axis).*100, flip(db(e_1)), 'LineWidth', 1.8);
xlabel('\delta [cm]');
ylabel('e_1 [dB]');
set(gca, 'FontSize', 12);

%% b) Evaluate the error e2 as a function of the length δ of the conical 
% sections and plot it in Matlab.

e_2 = zeros(1, length(n_sections));
[Z2_max, Z2_freqs] = findpeaks(abs(Z_2), f_axis, 'NPeaks', 5);
for i = 1:length(e_2)
    [Z1_max, Z1_freqs] = findpeaks(abs(Z_1(i, :)), f_axis, 'NPeaks', 5);
    e_2(1, i) = sum(abs(Z1_freqs-Z2_freqs));
end
figure();
stem(flip(delta_axis).*100, flip(e_2), 'LineWidth', 1.8);
xlabel('\delta [cm]');
ylabel('e_2');
set(gca, 'FontSize', 12);

[e2_min, n_sec_opt] = min(e_2);
delta_opt = l_exp_horn./n_sec_opt;
Z_1_opt = Z_1(n_sec_opt, :);

%% c) Compute the impedance of the approximated model when the radiation 
% load is kept into account and evaluate the error brought by neglecting 
% it in terms of the metric e1.

Z_l0 = 0.25.*w.^2.*(rho/(pi*c)) + 0.61*1i.*w*(rho/(pi*a));
r_1 = a_0*exp(m*(l_exp_horn-delta_opt));
x_1 = delta_opt*(r_1/(a-r_1));
flare_angle = atan(a/x_1);
S_s = (2*S_m)/(1+cos(flare_angle));
Z_l = (S_m/S_s).*Z_l0;
Z_1_opt_rad = approximate_exponential_horn(l_exp_horn, a_0, m, ...
    delta_opt, k, Z_l, rho, c);

figure();
plot(f_axis, db(abs(Z_1_opt)), 'LineWidth', 1.8);
hold on;
plot(f_axis, db(abs(Z_1_opt_rad)), 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('|Z| [dB]');
legend('Z_1', 'Z_1 with radiation load');
set(gca, 'FontSize', 12);

e_1c = mean(abs(Z_1_opt_rad-Z_2).^2, 'omitnan');
e_1c_db = db(e_1c);

%% d) Compute the impedance of the compound horn and list in a table the 
% frequencies of the first ten maxima of the impedance.

Z_comp_n = Z_1_opt_rad.*cos(k.*l_pipe)+1i.*Z_0.*sin(k.*l_pipe);
Z_comp_d = 1i.*Z_1_opt_rad.*sin(k.*l_pipe)+Z_0.*cos(k.*l_pipe);
Z_comp = Z_0*(Z_comp_n./Z_comp_d);

figure();
plot(f_axis, db(abs(Z_comp)), 'LineWidth', 1.8);
xlabel('f [Hz]');
ylabel('|Z_c| [dB]');
set(gca, 'FontSize', 12);

[Z_comp_max, Z_comp_freqs] = findpeaks(abs(Z_comp), f_axis, 'NPeaks', 10);
fmt = ['Maxima of the compund horn impedance:\n' repmat(' %1.3f\n', 1, ...
    numel(Z_comp_freqs))];
fprintf(fmt, Z_comp_freqs);

%% e) Evaluate the inharmonicity of the succession of impedance maxima for 
% the different setups: exponential horn, exponential horn with impedance 
% radiation and compound horn with impedance radiation.

[Z_1_opt_max, Z_1_opt_freqs] = findpeaks(abs(Z_1_opt), f_axis, 'NPeaks', 10);
fmt = ['Maxima of the approximated exponential horn impedance:\n' ...
    repmat(' %1.3f\n', 1, numel(Z_1_opt_freqs))];
fprintf(fmt, Z_1_opt_freqs);

[Z_1_opt_rad_max, Z_1_opt_rad_freqs] = findpeaks(abs(Z_1_opt_rad), f_axis, 'NPeaks', 10);
fmt = ['Maxima of the approximated exponential horn impedance with radiation load:\n' ...
    repmat(' %1.3f\n', 1, numel(Z_1_opt_rad_freqs))];
fprintf(fmt, Z_1_opt_rad_freqs);
