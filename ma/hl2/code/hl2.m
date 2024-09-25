clc;
clearvars;
close all;

%% Import paths
addpath('./libraries/polimi');
addpath('./functions');

%% Resonator parameters
V_0 = 0.1; % m^3
l = 0.1; % 10 cm
S = 100; % m^2
c = 343; % m/s
rho = 1.2; % Kg/m^3

%% Electric analogs
a = sqrt(S/pi);
end_correction_l = (0.61 + 8/(3*pi))*a;
l_tot = l + end_correction_l;
R = rho*c/S;
L = rho*l_tot/S;
C = V_0/(rho*c^2);

%% Simulation parameters
Fs = 44100;
Ts = 1/Fs;
sim_duration = '3.0';
fig_index = 1;

%% Ex.1: Model the response of a single Helmholtz resonator.

% (a) Set a simulation in Simscape and plot the frequency response
single_res_model = generate_helmholtz_tree_model(0, 1, 'single', Fs);
single_sim_out = sim(single_res_model, 'StopTime', sim_duration);
L_data = length(single_sim_out.flow_0_0.Data);
N = 2^nextpow2(L_data);
pressure = fft(squeeze(single_sim_out.pressure_0.Data), N)/L_data;
acoustic_flow = fft(squeeze(single_sim_out.flow_0_0.Data), N)/L_data;
single_res_frf = acoustic_flow./pressure;
f_axis = Fs/2*linspace(0, 1, N/2+1);
single_res_frf = single_res_frf(1:N/2+1);
single_res_frf(2:end) = 2*single_res_frf(2:end);
single_res_fig = figure(fig_index);
fig_index = fig_index+1;
ax_1 = axes('Parent', single_res_fig);
hold on;
plot(ax_1, f_axis, db(abs(single_res_frf)), 'LineWidth', 2);
hold off;
grid on;
xlim(ax_1, [0, 1000]);
ylim(ax_1, [-60, -5]);
xlabel(ax_1, 'f [Hz]', 'FontSize', 12);
ylabel(ax_1, '|H(\omega)| [dB]', 'FontSize', 12);

% (b) Compute the natural frequency of the resonator analytically and 
% verify that it matches the results of the simulation
damping_factor = R/2*sqrt(C/L);
q_factor = 1/R*sqrt(L/C);
w_0 = 1/sqrt(L*C);
w_d = w_0*sqrt(1-damping_factor^2);
res_freq_an = w_d/(2*pi);
[peaks, locs] = findpeaks(abs(single_res_frf));
res_freq_sim = f_axis(locs);
fprintf('Simulated natural frequency is %4.8f\n', res_freq_sim);
fprintf('Analytical natural frequency is %4.8f\n', res_freq_an);

%% Ex.2: Combine more resonators in a tree and analyze the response 
% obtained.

% (a) Use the RLC circuit defined in Ex.1 and connect its replicas to 
% build a NxK tree. Use the same parameters for each component and analyze
% the frequency response using as output the current in one of the leaves 
% (𝑈𝑛,𝑘) and pressure 𝑝0 as input


% (b) Try with different N and K and highlight what these parameters 
% control in the final response

K_pairs = [1 3; 2 3; 3 3; 4 3];
N_pairs = [4 1; 4 2; 4 3; 4 4];

%% Unbalanced tree simulation (K fixed)
unbal_k_fixed_fig = figure(fig_index);
fig_index = fig_index+1;
ax_2 = axes('Parent', unbal_k_fixed_fig);
plot(ax_2, f_axis, db(abs(single_res_frf)), 'LineWidth', 2);
hold on;
for i = 1:length(K_pairs)
    N_K_pair = K_pairs(i, :);
    sim_model = generate_helmholtz_tree_model( ...
        N_K_pair(1), N_K_pair(2), 'unbalanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    leaf_id = sprintf('_%d_%d', N_K_pair(1), 0);
    pressure_var = strcat('pressure', '_0');
    flow_var = strcat('flow', leaf_id);
    L_data = length(sim_out.(flow_var).Data);
    N = 2^nextpow2(L_data);
    pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
    acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
    freq_response = acoustic_flow./pressure;
    f_axis = Fs/2*linspace(0, 1, N/2+1);
    freq_response = freq_response(1:N/2+1);
    freq_response(2:end) = 2*freq_response(2:end);
    plot(ax_2, f_axis, db(abs(freq_response)), 'LineWidth', 2);
end
hold off;
grid on;
xlim(ax_2, [0, 1800]);
ylim(ax_2, [-70, -5]);
xlabel(ax_2, 'f [Hz]', 'FontSize', 12);
ylabel(ax_2, '|H(\omega)| [dB]', 'FontSize', 12);
legend(ax_2, '0x1','1x3','2x3','3x3', '4x3');

%% Unbalanced tree simulation (N fixed)
unbal_n_fixed_fig = figure(fig_index);
fig_index = fig_index+1;
ax_3 = axes('Parent', unbal_n_fixed_fig);
plot(ax_3, f_axis, db(abs(single_res_frf)), 'LineWidth', 2);
hold on;
for i = 1:length(N_pairs)
    N_K_pair = N_pairs(i, :);
    sim_model = generate_helmholtz_tree_model( ...
        N_K_pair(1), N_K_pair(2), 'unbalanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    leaf_id = sprintf('_%d_%d', N_K_pair(1)-1, N_K_pair(2)-1);
    pressure_var = strcat('pressure', '_0');
    flow_var = strcat('flow', leaf_id);
    L_data = length(sim_out.(flow_var).Data);
    N = 2^nextpow2(L_data);
    pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
    acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
    freq_response = acoustic_flow./pressure;
    f_axis = Fs/2*linspace(0, 1, N/2+1);
    freq_response = freq_response(1:N/2+1);
    freq_response(2:end) = 2*freq_response(2:end);
    plot(ax_3, f_axis, db(abs(freq_response)), 'LineWidth', 2);
end
hold off;
grid on;
xlim(ax_3, [0, 1800]);
ylim(ax_3, [-70, -5]);
xlabel(ax_3, 'f [Hz]', 'FontSize', 12);
ylabel(ax_3, '|H(\omega)| [dB]', 'FontSize', 12);
legend(ax_3, '0x1', '4x1', '4x2', '4x3', '4x4');

%% Unbalanced tree simulation (K fixed)
bal_k_fixed_fig = figure(fig_index);
fig_index = fig_index+1;
ax_4 = axes('Parent', bal_k_fixed_fig);
plot(ax_4, f_axis, db(abs(single_res_frf)), 'LineWidth', 2);
hold on;
for i = 1:length(K_pairs)
    N_K_pair = K_pairs(i, :);
    sim_model = generate_helmholtz_tree_model( ...
        N_K_pair(1), N_K_pair(2), 'unbalanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    leaf_id = sprintf('_%d_%d', N_K_pair(1), 0);
    pressure_var = strcat('pressure', '_0');
    flow_var = strcat('flow', leaf_id);
    L_data = length(sim_out.(flow_var).Data);
    N = 2^nextpow2(L_data);
    pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
    acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
    freq_response = acoustic_flow./pressure;
    f_axis = Fs/2*linspace(0, 1, N/2+1);
    freq_response = freq_response(1:N/2+1);
    freq_response(2:end) = 2*freq_response(2:end);
    plot(ax_4, f_axis, db(abs(freq_response)), 'LineWidth', 2);
end
hold off;
grid on;
xlim(ax_4, [0, 1800]);
ylim(ax_4, [-70, -5]);
xlabel(ax_4, 'f [Hz]', 'FontSize', 12);
ylabel(ax_4, '|H(\omega)| [dB]', 'FontSize', 12);
legend(ax_4, '0x1','1x3','2x3','3x3', '4x3');

%% Balanced tree simulation (N fixed)
bal_n_fixed_fig = figure(fig_index);
fig_index = fig_index+1;
ax_5 = axes('Parent', bal_n_fixed_fig);
plot(ax_5, f_axis, db(abs(single_res_frf)), 'LineWidth', 2);
hold on;
for i = 1:length(N_pairs)
    N_K_pair = N_pairs(i, :);
    sim_model = generate_helmholtz_tree_model( ...
        N_K_pair(1), N_K_pair(2), 'balanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    leaf_id = sprintf('_%d_%d', N_K_pair(1)-1, N_K_pair(2)-1);
    pressure_var = strcat('pressure', '_0');
    flow_var = strcat('flow', leaf_id);
    L_data = length(sim_out.(flow_var).Data);
    N = 2^nextpow2(L_data);
    pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
    acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
    freq_response = acoustic_flow./pressure;
    f_axis = Fs/2*linspace(0, 1, N/2+1);
    freq_response = freq_response(1:N/2+1);
    freq_response(2:end) = 2*freq_response(2:end);
    plot(ax_5, f_axis, db(abs(freq_response)), 'LineWidth', 2);
end
hold off;
grid on;
xlim(ax_5, [0, 1800]);
ylim(ax_5, [-70, -5]);
xlabel(ax_5, 'f [Hz]', 'FontSize', 12);
ylabel(ax_5, '|H(\omega)| [dB]', 'FontSize', 12);
legend(ax_5, '0x1', '4x1', '4x2', '4x3', '4x4');

%% (c) What happens if you change the location in which you evaluate the 
% frequency response inside the tree hierarchy? Elaborate on that, showing 
% some examples

level_pairs = [1 3; 2 3; 3 3; 4 3; 4 1; 4 2; 4 3; 4 4];

%% Balanced tree per-level simulation
level_fig = figure(fig_index);
set(0, 'CurrentFigure', level_fig);
fig_index = fig_index+1;
plot_pos = 1;
for i = 1:length(level_pairs)
    N_K_pair = level_pairs(i, :);
    sim_model = generate_helmholtz_tree_model(N_K_pair(1), N_K_pair(2), 'balanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    if (plot_pos == 5)
        level_fig = figure(fig_index);
        set(0, 'CurrentFigure', level_fig);
        fig_index = fig_index+1;
        plot_pos = 1;
    end
    ax_fig = subplot(2, 2, plot_pos);
    plot_pos = plot_pos + 1;
    for j = 0:N_K_pair(1)
        child_id = sprintf('_%d_%d', j, 0);
        pressure_var = strcat('pressure', '_0');
        flow_var = strcat('flow', child_id);
        L_data = length(sim_out.(flow_var).Data);
        N = 2^nextpow2(L_data);
        pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
        acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
        freq_response = acoustic_flow./pressure;
        f_axis = Fs/2*linspace(0, 1, N/2+1);
        freq_response = freq_response(1:N/2+1);
        freq_response(2:end) = 2*freq_response(2:end);
        plot(ax_fig, f_axis, db(abs(freq_response)), 'LineWidth', 2);
        hold on;
    end
    hold off;
    grid on;
    xlim(ax_fig, [0, 1800]);
    ylim(ax_fig, [-100, -5]);
    xlabel(ax_fig, 'f [Hz]', 'FontSize', 12);
    ylabel(ax_fig, '|H(\omega)| [dB]', 'FontSize', 12);
    title(sprintf('(%d, %d)',  N_K_pair(1),  N_K_pair(2)));
end

%% Unbalanced tree per-level simulation
level_fig = figure(fig_index);
set(0, 'CurrentFigure', level_fig);
fig_index = fig_index+1;
plot_pos = 1;
for i = 1:length(level_pairs)
    N_K_pair = level_pairs(i, :);
    sim_model = generate_helmholtz_tree_model(N_K_pair(1), N_K_pair(2), 'unbalanced', Fs);
    sim_out = sim(sim_model, 'StopTime', sim_duration);
    if (plot_pos == 5)
        level_fig = figure(fig_index);
        set(0, 'CurrentFigure', level_fig);
        fig_index = fig_index+1;
        plot_pos = 1;
    end
    ax_fig = subplot(2, 2, plot_pos);
    plot_pos = plot_pos + 1;
    for j = 0:N_K_pair(1)
        child_id = sprintf('_%d_%d', j, 0);
        pressure_var = strcat('pressure', '_0');
        flow_var = strcat('flow', child_id);
        L_data = length(sim_out.(flow_var).Data);
        N = 2^nextpow2(L_data);
        pressure = fft(squeeze(sim_out.(pressure_var).Data), N)/L_data;
        acoustic_flow = fft(squeeze(sim_out.(flow_var).Data), N)/L_data;
        freq_response = acoustic_flow./pressure;
        f_axis = Fs/2*linspace(0, 1, N/2+1);
        freq_response = freq_response(1:N/2+1);
        freq_response(2:end) = 2*freq_response(2:end);
        plot(ax_fig, f_axis, db(abs(freq_response)), 'LineWidth', 2);
        hold on;
    end
    hold off;
    grid on;
    xlim(ax_fig, [0, 1800]);
    ylim(ax_fig, [-100, -5]);
    xlabel(ax_fig, 'f [Hz]', 'FontSize', 12);
    ylabel(ax_fig, '|H(\omega)| [dB]', 'FontSize', 12);
    title(sprintf('(%d, %d)',  N_K_pair(1),  N_K_pair(2)));
end
