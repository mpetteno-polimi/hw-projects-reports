clc;
close all;
clearvars;
addpath(genpath('functions'));

%% Parameters
% Parameters circular membrane
r_m = 0.15;           % Radius = 0.15 m
T_m = 10;             % Tension = 10 N/m
sigma_m = 0.07;       % Unit surface weight σ=0.07 kg/m^2
Q_m = 25;

% Parameters circular plate
r_p = r_m;            % Radius = 0.15 m
h_p = 0.001;          % Thickness = 1mm
E_p = 69e9;           % Young's modulus = 69GPa
rho_p = 2700;         % Volume density = 2700kg/m^3
nu_p = 0.334;         % Poisson's ratio = 0.334
Q_p = 50; 

% Parameters iron string
rho_s = 5000;         % Volume density of the string rho = 5000 kg/m^3
r_s = 0.001;          % Radius of the string = 0.001 m
l_s = 0.4;            % Length of the string = 0.4 m
E_s = 200e9;          % Young's modulus = 200GPa


%% PART 1 - Circular membrane characterization

%% a) Propagation speed in the membrane
c_m = sqrt(T_m/sigma_m); % [m/s]

%% b) Frequency of the first eighteen modes for this membrane (ordered by frequency)
n_modes_m = 18;   

% Bessel function zeros
n = (0:9)';
k = 10;
kind = 1;
besselz = besselzero(n, k, kind);
sorted_bessel_zeros = sort(besselz(:));
bessel_zeros = sorted_bessel_zeros(1:n_modes_m);
bessel_zeros_idx = zeros(n_modes_m, 2);
for i=1:n_modes_m
    [m, n] = find(besselz == bessel_zeros(i));
    bessel_zeros_idx(i, :) = [m - 1, n];
end

% Frequencies
mode_k_m = bessel_zeros/r_m;
mode_w_m = mode_k_m*c_m;
mode_f_m = mode_w_m/(2*pi);

% Modeshapes
A = 1;
radial_res = 100;
angular_res = 100;
phi = linspace(0, 2*pi, angular_res);
a = linspace(0, r_m, radial_res);
[R, PHI] = meshgrid(a, phi);
mode_shapes_m = zeros(radial_res, angular_res, n_modes_m);
for i=1:n_modes_m
    idx = bessel_zeros_idx(i, :);
    m = idx(1);
    k = mode_k_m(i);
    mode_shapes_m(:, :, i) = A.*exp(1i.*m.*PHI).*besselj(m, k.*R);
end

modes_to_plot = 6;
[x, y] = pol2cart(PHI, R);
figure();
for i=1:modes_to_plot
    subplot(modes_to_plot/3, 3, i);
    idx = bessel_zeros_idx(i, :);
    m = idx(1);
    n = idx(2);
    surf(x, y, real(mode_shapes_m(:, :, i)));
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title(sprintf('(%d, %d)', m, n));
    set(gca, 'FontSize', 12);
end

%% c) FRF modal superposition analysis
Fs = 48000;
T = 1/Fs;
duration = 2;
N = Fs*duration;
f_axis = 0:Fs/N:Fs-Fs/N;
w_axis = 2*pi*f_axis;
ex_r = 0.075;
ex_phi = deg2rad(0);
meas_r = 0.075;
meas_phi = deg2rad(195);
damping_loss_m = 1/Q_m;

% compute receptance
receptance_m = zeros(1, N);
for i = 1:n_modes_m
    fprintf('----------------- MODE %d -----------------------\n', i);
    idx = bessel_zeros_idx(i, :);
    m = idx(1);
    n = idx(2);
    k = mode_k_m(i);
    w = mode_w_m(i);
    fprintf('circular nodes = %d\n', n);
    fprintf('diametrical nodes = %d\n', m);
    fprintf('frequency = %4.8f\n', w/(2*pi));
    % compute modal mass for normalization
    polarfun = @(phi,r) r*sigma_m.*abs((A.*exp(1i.*m.*phi).*besselj(m, k.*r))).^2;
    modal_mass = integral2(polarfun, 0, 2*pi, 0, r_m);
    % compute excitation and measurement points displacements
    Z_ex = A.*exp(1i.*m.*ex_phi).*besselj(m, k.*ex_r);
    fprintf('exitaction = %4.8f\n', Z_ex);
    Z_meas = A.*exp(1i.*m.*meas_phi).*besselj(m, k.*meas_r);
    fprintf('measurement = %4.8f\n', Z_ex);
    % compute the receptance for the mode
    mode_receptance = Z_ex*Z_meas./(modal_mass*(-w_axis.^2+w^2*(1+1i*damping_loss_m)));
    receptance_m = receptance_m + mode_receptance;
end

% input force plot
t_axis = (0:N-1)*T;
input_force = 0.1.*exp(-(t_axis - 0.03).^2./(0.01)^2);
figure();
plot(t_axis, input_force, 'LineWidth', 1.8, 'Color', '#0072BD');
xlabel('Time (s)');
ylabel('Force (N)');
set(gca, 'FontSize', 20);

% receptance plots
figure();
plot(f_axis, 20*log10(abs(receptance_m)), 'LineWidth', 1.8, 'Color', '#0072BD');
xlabel('Frequency (Hz)');
ylabel('Magnitude (db)');
set(gca, 'FontSize', 20);

ifft_receptance_m = ifft(receptance_m, 'symmetric');
figure();
plot(t_axis, ifft_receptance_m, 'LineWidth', 1.8, 'Color', '#0072BD');
xlabel('Time (s)');
ylabel('m/N');
set(gca, 'FontSize', 20);

% sensor displacement plot
sensor_displacement = conv(input_force, ifft_receptance_m);
displacement_dur = 2*Fs;
conv_axis = (0:displacement_dur-1)/Fs;
figure();
plot(conv_axis, sensor_displacement(1:displacement_dur), 'LineWidth', 1.8, 'Color', '#0072BD');
xlabel('Time (s)');
ylabel('Displacement (m)');
set(gca, 'FontSize', 20);

%% PART 2 - Circular plate (alluminium) characterization

%% d) Propagation speed of quasi-longitudinal and longitudinal waves
c_qlw = sqrt(E_p/(rho_p*(1-(nu_p^2)))); % quasi-longitudinal waves spped
c_lw = sqrt(E_p*(1-nu_p)/(rho_p*(1+nu_p)*(1-2*nu_p))); % longitudinal waves speed

%% e) Plot of propagation speed of bending waves as function of the frequency
freq_axis = linspace(0, 10e3, 10e6);
omega = freq_axis.*(2*pi);
k = sqrt((omega.*sqrt(12))/(c_qlw*h_p));
v = omega./k; % propagation speed of bending waves

figure();
plot(freq_axis, v, 'LineWidth', 1.8, 'Color', '#0072BD');
hold on;
grid on;
xlabel('Frequency [Hz]');
ylabel('Propagation speed [m/s]');
set(gca, 'FontSize', 20);

%% f) Modal frequencies of the first 5 bending modes of the plate
n_modes_p = 5;
mode_coeff_p = [1; 2.08; 3.41; 3.89; 5.00];
mode_f_p = zeros(1, n_modes_p);
f_01_p = (0.4694*c_qlw*h_p)/(r_p^2);
for n = 1:n_modes_p
    mode_f_p(1, n) = mode_coeff_p(n)*f_01_p;
end


%% PART 3 - Interaction between coupled systems

%% g) Tension of the string
% c = sqrt(T/sigma) -> sigma is linear density
% n-th mode of the string -> w_n = n*c*pi/l_s -> foundamental mode -> n = 1
% w_01_p = c*pi/l_s = sqrt(T/sigma)*pi/l_s -> T = (w_01_p*l_s/pi)^2*sigma_s
w_01_p = f_01_p*2*pi;
sigma_s = pi*r_s^2*rho_s;
T_s = (w_01_p*l_s/pi)^2*sigma_s;

%% h) Frequencies of the first five modes of the string considering stiffness
n_modes_s = 5;
K = r_s/2;
S = pi*r_s^2;
B = pi*E_s*S*K^2/(T_s*l_s^2);
mode_f_stiff_s = zeros(1, n_modes_s);
clamped_edge_factor = 1+2*sqrt(B)/pi+4*B/pi^2;
for n = 1:n_modes_s
    mode_f_stiff_s(1, n) = n*f_01_p*sqrt(1+B*n^2) * clamped_edge_factor;
end

%% i) Frequencies of the modes of the string-soundboard system considering coupliing between the plate and the string
n_modes_coupl = 5;
mode_f_s = f_01_p.*linspace(1, 5, n_modes_coupl); % string natural modes without stiffness
m_s = sigma_s*l_s;
m_p = rho_p*h_p*pi*r_p^2;
weak_strong_coupl_bound = pi^2/(4*Q_p^2);
fprintf('Weak-Strong coefficient threshold %4.8f\n', weak_strong_coupl_bound);
for n = 1:n_modes_coupl
    fprintf('----------------- MODE %d -----------------------\n', n);
    mode_coupling_coeff = m_s/(n^2*m_p);
    fprintf('Coupling coefficient %4.8f\n', mode_coupling_coeff);
    f_s = mode_f_s(1, n);
    f_p = mode_f_p(1, n);
    if (mode_coupling_coeff < weak_strong_coupl_bound)
        fprintf('Weak coupling\n');
        fprintf('\tString mode is %4.2f Hz, plate mode is %4.2f Hz\n', f_s, f_p);
    else
        fprintf('Strong coupling\n');
        fprintf('\tString mode is %4.2f Hz, plate mode is %4.2f Hz\n', f_s, f_p);
        if (f_s == f_p)
            fprintf('GRAFICO 1:\n Q factor = %4.0f\n Coupling coefficient (asse x) = %4.8f\n' , Q_p, mode_coupling_coeff);
            prompt = 'Enter value of (\Omega_+-\Omega_-)/\omega: ';
            Strong_equal_y = input(prompt);
            delta = Strong_equal_y*f_p/2;
            Freq_plus = f_p + delta;
            Freq_minus = f_p - delta;
            fprintf('Resonance frequencies: %4.2f Hz and %4.2f Hz\n', Freq_plus, Freq_minus);
        else
            Strong_diff_x = (f_s - f_p)/f_p;
            fprintf('GRAFICO 2 (Strong): asse x -> %4.8f\n', Strong_diff_x);
            prompt1 = 'Enter first value of (\Omega_(1)-\omega_p)/\omega_p: ';
            Strong_diff_y_1 = input(prompt1);
            Freq_uno = (Strong_diff_y_1+1)*f_p;
            prompt2 = 'Enter second value of (\Omega_(2)-\omega_p)/\omega_p: ';
            Strong_diff_y_2 = input(prompt2);
            Freq_due = (Strong_diff_y_2+1)*f_p;
            fprintf('Resonance frequencies: %4.2f Hz and %4.2f Hz\n', Freq_uno, Freq_due);
        end
    end
end