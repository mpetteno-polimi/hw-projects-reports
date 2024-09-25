%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6
% Complete model of a guitar
% This script call the simulink file that implements the simple electrical
% equivalent of a guitar. The simulation is executed and the resulting
% current is resampled with the desired sampling frequency. Finally the
% sound is plotted in time and frequency domain and saved on disk. 
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

addpath('./functions');
simulink_folder = './simulink';         % Simulink projects folder
addpath(simulink_folder);

%% Setup

fs = 44100;                             % Sampling frequency
signalLen = 5;                          % Signal length
t = 0:1/fs:signalLen-1/fs;              % Time axis

%% Simulation (Without string model)
% run the simulink simulation using the command sim (see doc sim).
sim('ex_6_volt_gen'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I_vg = resample(I1, t);

%% Simulation (With string model)
% run the simulink simulation using the command sim (see doc sim).
sim('ex_6_string_model'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I_sm = resample(I2, t);

%% Plot
% Plot the resampled signals in time
fig_1 = figure(1);
subplot(2, 1, 1);
plot(I_vg.time, I_vg.data);
xlabel('Time [s]');
ylabel('Voltage generator');
subplot(2, 1, 2);
plot(I_sm.time, I_sm.data);
xlabel('Time [s]');
ylabel('String model');

% Normalize the signals
soundWave_vg = I_vg.data;
soundWave_vg = soundWave_vg./max(abs(soundWave_vg));
soundWave_sm = I_sm.data;
soundWave_sm = soundWave_sm./max(abs(soundWave_sm));

% Plot the signals frequency content as magnitude and phase
% Voltage generator model
fig_2 = figure(2);
[S_vg, magS_vg, angleS_vg, f_vg, df_vg] = myFFT(soundWave_vg, fs);
plotFFT_linearFreqScale(magS_vg, angleS_vg, f_vg, df_vg, fs, 2000, fig_2);
% String model
fig_3 = figure(3);
[S_sm, magS_sm, angleS_sm, f_sm, df_sm] = myFFT(soundWave_sm, fs);
plotFFT_linearFreqScale(magS_sm, angleS_sm, f_sm, df_sm, fs, 2000, fig_3);

fig_4 = figure(4);
subplot(2, 1, 1);
plot(t, soundWave_vg);
xlabel('Time [s]');
ylabel('Voltage generator');
subplot(2, 1, 2);
plot(t, soundWave_sm);
xlabel('Time [s]');
ylabel('String model');

%% Play
% Play the sound
sound(soundWave_vg, fs);
pause(3);
sound(soundWave_sm, fs);
 
%% Save
% Save on disk
disp('Save file on disk...')
id_number = '10868930';
surname = 'petteno';
fileName = sprintf('%s_%s_complete_guitar_vg.wav', id_number, surname);
audiowrite(fileName, soundWave_vg, fs);
fileName = sprintf('%s_%s_complete_guitar_sm.wav', id_number, surname);
audiowrite(fileName, soundWave_sm, fs);
