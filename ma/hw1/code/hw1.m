clc;
close all;
clearvars;

%% Parameters

% Internal volume
resonator_heigth = 0.25; % 25cm
resonator_width = 0.25; % 25cm
resonator_length = 0.18; % 18cm
resonator_volume = resonator_heigth * resonator_width * resonator_length;

% Resonator's neck
neck_length = 0.06; % 0.6cm
neck_radius = 0.025; % 2.5cm
neck_cross_section = pi*neck_radius^2;

% Air density
absolute_pressure = 101.325e3; % at sea level 1 atm = 101.325 kPa
absolute_temperature = 293.15; % 20 °C = 293.15K
dry_air_molecular_mass = 4.81e-26;
k = physconst('Boltzmann');
air_density = (absolute_pressure*dry_air_molecular_mass)/(k*absolute_temperature);

% Speed of sound
air_bulk_modulus = 1.4e5; % 1.4*10^5 Pa
c = sqrt(air_bulk_modulus/air_density);

%% 1) Resonance frequency (without virtual neck elongation)

neck_air_mass_no_vn = air_density*neck_cross_section*neck_length;
spring_constant_no_vn = air_density*(neck_cross_section*c)^2/resonator_volume;
w0_no_vn = sqrt(spring_constant_no_vn/neck_air_mass_no_vn);
f0_no_vn = w0_no_vn/(2*pi);
lambda_no_vn = c/f0_no_vn;

%% 2) Resonance frequency (with virtual neck elongation)

vn_end_correction = neck_radius*0.61; % unflanged (one side)
neck_length_vn = neck_length + vn_end_correction;
neck_air_mass_vn = air_density*neck_cross_section*neck_length_vn;
spring_constant_vn = spring_constant_no_vn;
w0_vn = sqrt(spring_constant_vn/neck_air_mass_vn);
f0_vn = w0_vn/(2*pi);
lambda_vn = c/f0_no_vn;

%% 3) Critically damped system's resistance
% a = R/2m -> if crit. damp. -> a = w0 -> R = 2m*w0;

R_no_vn = 2*neck_air_mass_no_vn*w0_no_vn;
R_vn = 2*neck_air_mass_vn*w0_vn;

%% 4) System with virtual neck impedance and R = 5*10^-4 kg/s
R = 5e-4;
freqs = linspace(81.5, 84);
w = 2*pi.*freqs;
X = w.*neck_air_mass_vn - spring_constant_vn./w;
Z = R + 1i*X;
figure();
real_plot = semilogx(freqs, real(Z), LineStyle='-', LineWidth=1.8, DisplayName='R');
hold on;
imag_plot = semilogx(freqs, imag(Z), LineStyle='-', LineWidth=1.8, DisplayName='X');
hold on;
abs_plot = semilogx(freqs, abs(Z), LineStyle='--', LineWidth=1.8, DisplayName='|Z|');
hold on;
xline(f0_vn, LineStyle='--', LineWidth=1.2, label={'Resonance', f0_vn}, ...
    LabelOrientation='horizontal');
hold on;
semilogx(freqs, zeros(length(freqs), 1), Color='k');
hold off;
xlabel('frequency');
ylabel('impedance');
legend([real_plot, imag_plot, abs_plot], Location='northeastoutside');
set(gca, 'box', 'off');

%% 5) Q factor and time decay
Q_factor = spring_constant_vn/(w0_vn*R);
decay_time = 2*neck_air_mass_vn/R;

%% 6) Q and resonance frequency as function of R = [0:0.5] kg/s
R_axis = linspace(0, 0.5, 1e3);
alpha = R_axis./(2*neck_air_mass_vn);
w0_func = sqrt(w0_vn^2 - alpha.^2);
f0_func = w0_func./(2*pi);
Q_factor_func = spring_constant_vn./(w0_func.*R_axis);
figure();
plot(R_axis, Q_factor_func, LineStyle='-', LineWidth=1.8, DisplayName='Q factor');
hold on;
plot(R_axis, w0_func, LineStyle='-', LineWidth=1.8, DisplayName='Damped resonance frequency');
hold on;
yline(w0_vn, Color='g', LineStyle='--', LineWidth=1.5, DisplayName='Natural resonance frequency');
hold on;
plot(R_axis, alpha, LineStyle=':', LineWidth=1.5, DisplayName='alpha');
hold off;
xlabel('resistance');
legend;
set(gca, 'box', 'off');
