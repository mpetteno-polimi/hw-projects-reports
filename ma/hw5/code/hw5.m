clear;
close all;
clc;

addpath("./functions");

%% Parameters

% Air density
abs_p = 101.325e3; % at sea level 1 atm = 101.325 kPa
abs_temp = 293.15; % 20 °C = 293.15K
dry_air_mol_mass = 4.81e-26;
k_boltz = physconst('Boltzmann');
rho = (abs_p*dry_air_mol_mass)/(k_boltz*abs_temp); % Kg/m^3

% Air kinematic viscosity
nu = 1.5e-5; % m^2/s

% Speed of sound
air_bulk_modulus = 1.4e5; % 1.4*10^5 Pa
c = sqrt(air_bulk_modulus/rho);

% Sampling
N = 2^16;
f_min = 0;
f_max = 2000;
f_axis = f_min:f_max/N:f_max-f_max/N;
w = 2*pi*f_axis;
k = w./c;

%% First component: resonator
% The resonator is shaped as a cone whose conical semiangle is 0.75°. The 
% instrument is aimed at being a treble recorder, with a length of 0.45m. 
% For the sake of simplicity, we consider that only two finger holes are 
% present.
theta = deg2rad(0.75);
l_r = 0.45;

%% Question 1 - Find the diameter of the cone at the resonator head and 
% foot so that the note produced when all the finger holes are closed is 
% E4 (329.63 Hz).
f_0 = 329.63;
w_0 = 2*pi*f_0;
k_0 = w_0/c;
% Bore radius boundaries
x_f_min = l_r;
x_f_max = 7*l_r/3;
r_f_min = x_f_min*tan(theta);
r_f_max = x_f_max*tan(theta);
% Expanded range bore impedance
x_f = linspace(0, 5, N);
r_f = x_f.*tan(theta);
r_h = (x_f+l_r).*tan(theta);
l_m = 0.6.*r_h; % mouth end correction
l_f = 0.6.*r_f; % foot end correction
L = l_r + l_m + l_f; % total length
Z = conical_horn_impedance(L, k_0, r_f, r_h, 0, rho, c);
Y = 1./Z;
% Plot bore admittaconical_horn_impedancence as function of x_f
figure();
plot(r_f.*100, db(abs(Y)), 'LineWidth', 1.8);
hold on;
xline(r_f_min*100, '--r', sprintf('r_{min} = %.2f cm', r_f_min*100), ...
    'LabelOrientation', 'aligned');
hold on;
xline(r_f_max*100, '--r', sprintf('r_{max} = %.2f cm', r_f_max*100), ...
    'LabelOrientation', 'aligned');
xlabel('r_f [cm]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);
% Valid range bore impedance
x_f = linspace(x_f_min, x_f_max, N);
r_f = x_f.*tan(theta);
r_h = (x_f+l_r).*tan(theta);
l_m = 0.6.*r_h; % mouth end correction
l_f = 0.6.*r_f; % foot end correction
L = l_r + l_m + l_f; % total length
Z = conical_horn_impedance(L, k_0, r_f, r_h, 0, rho, c);
Y = 1./Z;
% Optimum bore parameters
[Y_max, Y_max_idx] = max(abs(Y));
x_f_opt = x_f(Y_max_idx);
r_f_opt = x_f_opt*tan(theta);
r_h_opt = (x_f_opt+l_r)*tan(theta);
l_m_opt = 0.6*r_h_opt; % mouth end correction
l_f_opt = 0.6*r_f_opt; % foot end correction
L_opt = l_r + l_m_opt + l_f_opt; % total length
Z_opt = conical_horn_impedance(L_opt, k, r_f_opt, r_h_opt, 0, rho, c);
Y_opt = 1./Z_opt;
% Plot optimum bore impedance
figure();
plot(f_axis, db(abs(Y_opt)), 'LineWidth', 1.8);
hold on;
xline(f_0, '--r', sprintf('E4 = %.2f Hz', f_0), ...
    'LabelOrientation', 'horizontal');
xlabel('f [Hz]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);

%% Question 2 - Find the position of the last finger hole (i.e. the one 
% closest to the resonator foot) in order to produce the note F4# 
% (369.99 Hz) when it is open. Consider the finger hole diameter to be 
% equal to the bore diameter at the resonator foot (simplification).
f_1 = 369.99;
w_1 = 2*pi*f_1;
k_1 = w_1/c;
% Admittance(D)
D_axis = linspace(0, l_r, N);
% Pipe 2
L_p2 = D_axis + l_f_opt;
r_m = (x_f_opt+D_axis).*tan(theta);
Z_p2 = conical_horn_impedance(L_p2, k_1, r_f_opt, r_m, 0, rho, c);
Z_h = cylindrical_pipe_impedance(l_f_opt, k_1, r_f_opt, rho, c);
Z_l = Z_p2.*Z_h./(Z_p2+Z_h);
% Pipe 1
L_p1 = l_r - D_axis + l_m_opt;
Z_p1 = conical_horn_impedance(L_p1, k_1, r_m, r_h_opt, Z_l, rho, c);
Z_1 = Z_p1 + Z_l;
Y_1 = 1./Z_1;
% Plot admittance as function of D
figure();
plot(D_axis.*100, db(abs(Y_1)), 'LineWidth', 1.8);
xlabel('D [cm]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);
% Optimum bore parameters
[Y_1_max, Y_1_max_idx] = max(abs(Y_1));
D_1_opt = D_axis(Y_1_max_idx);
L_p2_opt = D_1_opt + l_f_opt;
r_m_opt = (x_f_opt+D_1_opt).*tan(theta);
Z_p2_opt = conical_horn_impedance(L_p2_opt, k, r_f_opt, r_m_opt, 0, rho, c);
Z_l_opt = Z_p2_opt.*Z_h./(Z_p2_opt+Z_h);
L_p1_opt = l_r - D_1_opt + l_m_opt;
Z_p1_opt = conical_horn_impedance(L_p1_opt, k, r_m_opt, r_h_opt, Z_l_opt, rho, c);
Z_1_opt = Z_p1_opt + Z_l_opt;
Y_1_opt = 1./Z_1_opt;
% Plot optimum bore impedance
figure();
plot(f_axis, db(abs(Y_1_opt)), 'LineWidth', 1.8);
hold on;
xline(f_1, '--r', sprintf('F#4 = %.2f Hz', f_1), ...
    'LabelOrientation', 'horizontal');
xlabel('f [Hz]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);

%% Question 3 - Find the position of the second last finger hole in order 
% to produce the note G4# (415.3 Hz) when the two finger holes are open. 
% Consider the finger hole diameter to be equal to the bore diameter at 
% the resonator foot (simplification).
f_2 = 415.30;
w_2 = 2*pi*f_2;
k_2 = w_2/c;
% Admittance(D)
D_axis = linspace(0, l_r, N);
% Pipe 3
L_p3 = D_1_opt + l_f_opt;
r_m_2 = (x_f_opt+D_1_opt).*tan(theta); 
Z_p3 = conical_horn_impedance(L_p3, k_2, r_m_2, r_f_opt, 0, rho, c);
% Pipe 2
L_p2 = D_axis - D_1_opt;
Z_l_1 = Z_p3.*Z_h./(Z_p3+Z_h);
r_m_1 = (x_f_opt+D_axis).*tan(theta);
Z_p2 = conical_horn_impedance(L_p2, k_2, r_m_1, r_m_2, Z_l_1, rho, c);
% Pipe 1
L_p1 = l_r - D_axis + l_m_opt;
Z_d = Z_p2 + Z_l_1;
Z_l_2 = Z_d.*Z_h./(Z_d+Z_h);
Z_p1 = conical_horn_impedance(L_p1, k_2, r_m_2, r_h_opt, Z_l_2, rho, c);
% Total impedance
Z_2 = Z_p1 + Z_l_2;
Y_2 = 1./Z_2;
% Plot admittance as function of D
figure();
plot(D_axis*100, db(abs(Y_2)), 'LineWidth', 1.8);
xlabel('D [cm]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);
% Optimum bore parameters
[Y_2_max, Y_2_max_idx] = max(abs(Y_2));
D_2_opt = D_axis(Y_2_max_idx);
% Pipe 3
L_p3_opt = D_1_opt + l_f_opt;
r_m_2_opt = (x_f_opt+D_1_opt).*tan(theta);
Z_p3_opt = conical_horn_impedance(L_p3_opt, k, r_m_2_opt, r_f_opt, 0, rho, c);
% Pipe 2
L_p2_opt = D_2_opt - D_1_opt;
Z_l_1_opt = Z_p3_opt.*Z_h./(Z_p3_opt+Z_h);
r_m_1_opt = (x_f_opt+D_2_opt).*tan(theta);
Z_p2_opt = conical_horn_impedance(L_p2_opt, k, r_m_1_opt, r_m_2_opt, Z_l_1_opt, rho, c);
% Pipe 1
L_p1_opt = l_r - D_2_opt + l_m_opt;
Z_d_opt = Z_p2_opt + Z_l_1_opt;
Z_l_2_opt = Z_d_opt.*Z_h./(Z_d_opt+Z_h);
Z_p1_opt = conical_horn_impedance(L_p1_opt, k, r_m_2_opt, r_h_opt, Z_l_2_opt, rho, c);
% Total impedance
Z_2_opt = Z_p1_opt + Z_l_2_opt;
Y_2_opt = 1./Z_2_opt;
% Plot optimum bore impedance
figure();
plot(f_axis, db(abs(Y_2_opt)), 'LineWidth', 1.8);
hold on;
xline(f_2, '--r', sprintf('G#4 = %.2f Hz', f_2), ...
    'LabelOrientation', 'horizontal');
xlabel('f [Hz]');
ylabel('|Y| [dB]');
set(gca, 'FontSize', 12);

%% Second component: flue channel
% The instrument is aimed at producing a spectrum whose centroid is at 1.7 
% kHz when the pressure difference between the player mouth and the flue 
% channel entrance is 55 Pa.
spectral_centroid = 1.7e3;
pressure_diff = 55;

% Question 4 - Find the flue channel thickness that complies with the 
% above specifications. For this pair of thickness and jet velocity, 
% compute the Reynolds number Re and assess the jet regime that is 
% undergoing at the flue channel exit and in the instrument mouth 
% (laminar, turbulent, etc).
U_j = sqrt(2*pressure_diff/rho); % Bernoulli approx. (constant jet vel.)
h = 0.3*U_j/spectral_centroid; % strongest amplification at f = 0.3*U_j/h
R_e = U_j*h/nu;
W = 4*h;

% Question 5 - Consider that the flue channel length is 25 mm. Find the 
% thickness of the boundary layer at the flue channel exit for the 
% specifications defined above (Question 4).
flue_channel_length = 0.025;
delta_l = sqrt(nu*flue_channel_length/U_j);

