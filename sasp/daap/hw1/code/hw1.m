clc; clear; close all;

%% Your name(s), student ID number(s)
%-------------------------------------------------------------------------%
% Matteo Pettenò, 10868930
% Marco Furio Colombo, 10537094
%-------------------------------------------------------------------------%

%% Here are some MATLAB function that you might find useful:
% audioread, soundsc, flipud, fliplr, xcorr, eig, eigs, filter, toeplitz,
% fft, ifft, pause, disp, ...

%% Read the wav files 
% load the source signal x
[x, Fx] = audioread('x.wav');

% load the microphone signal y
[y, Fy] = audioread('y.wav');

%% Parameters
% number of filter taps
M = 4000;

% length of the source signal x
N = size(x, 1);

%% Wiener-Hopf solution

% compute the autocorrelation matrix
maxlag = M-1;
[Rxx, Rxx_lags] = xcorr(x, x, maxlag, 'biased');
Rxx = Rxx(Rxx_lags >= 0);
R = toeplitz(Rxx);

% compute the cross-correlation vector
x_cross_padded = padarray(x, [length(y) - N, 0], 'post');
[Rxy, Rxy_lags] = xcorr(x_cross_padded, y, maxlag, 'biased');
Rxy = Rxy(Rxy_lags <= 0);
p = flipud(Rxy);

% compute the Wiener-Hopf solution
w_o = R\p;

%% Steepest Descent

% determine the step-size parameter
lambda_max = eigs(R, 1, 'largestabs');
mu_factor = 0.95;
mu = mu_factor * 2 / lambda_max;

% determine the global time constant of the Steepest Descent algorithm
lambda_min = eigs(R, 1, 'smallestabs');
tau = 1 / (2*mu*lambda_min);

% initialize the vector of filter taps
w = zeros(M, 1);

% perform the iterative update
iterations = 2000;
for i = 1:iterations
    w = w + mu*(p - R*w);
    
    % Optional: try with step-size optimization. 
    % With the same number of iterations, step-size optimization leads to 
    % better results but computational cost increases significantly.
    % mu = (p-R*w)'*(p-R*w)/((p-R*w)'*R*(p-R*w));

end

% compute the theoretical minimum MSE
J_min = var(y) - p'*w_o;

% compute the MSE associated with the current taps estimate
J_w = J_min + (w - w_o)'*R*(w - w_o);

% compute and discuss the cost functions ratio 
ratio = J_w / J_min;

%% Apply time-domain filter
y_hat_td = filter(w, 1, x);
%y_hat_td_opt = filter(w_o, 1, x); % for comparison

%% playback: source signal x
soundsc(x, Fx);
pause(2);

%% playback: microphone signal y
soundsc(y, Fy);
pause(2);

%% playback: time-domain filtered signal y_hat
soundsc(y_hat_td, Fx);
pause(2);
% for comparison
% playback: time-domain filtered signal y_hat_opt
%soundsc(y_hat_td_opt, Fx);
%pause(2);

%% Filtering in the frequency domain %
% determine the length of the signal after zero-padding
y_length = N + M - 1;
L = 2^nextpow2(y_length);

% compute the spectra
w_padded = padarray(w, [L - M, 0], 'post');
x_padded = padarray(x, [L - N, 0], 'post');
W = fft(w_padded);
X = fft(x_padded);

% perform frequency-domain filtering
Y = X.*W;

% transform back to time domain
y_hat_fd = ifft(Y);

%% playback: frequncy-domain filtered signal
soundsc(y_hat_fd, Fx);
pause(2);

%% OLA Filtering
% window length
wlen = 256;

% hop-size (must respect the COLA condition)
hop = wlen / 2; % COLA-C - could also be w_len/4

% define a tapered window
win = hann(wlen);

% determine the length of the windowed signal after zero-padding
L_ola = M + wlen - 1;

% compute the fft of the Wiener filter
W = fft(w, L_ola);

% compute the total number of frames to be processed via OLA
n_frames = ceil((N-wlen)/hop + 1);
% pad x to match the total numbers of frame
last_frame_end = (n_frames+1)*hop;
x = padarray(x, [last_frame_end-N, 0], 'post');

% initialize the output filtered signal
y_hat_ola = zeros(N + L_ola, 1);

% implement the OLA algorithm
for m = 0:(n_frames-1)
    shift = m*hop;
    x_m_frame_indexes = (shift+1):(shift+wlen); % indices for the mth frame
    x_m = x(x_m_frame_indexes) .* win;  % windowed mth frame
    x_m_pad_size = L_ola-wlen;
    x_m_padded = padarray(x_m, [x_m_pad_size, 0], 'post'); % zero pad the frame
    X_m = fft(x_m_padded); % mth frame FFT
    Y_m = X_m .* W; % frequency domain convolution
    y_m = ifft(Y_m); % filtered mth frame
    y_m_frame_indexes = (shift+1):(shift+L_ola);
    y_hat_ola(y_m_frame_indexes) = y_hat_ola(y_m_frame_indexes) + y_m; % overlap and add
end
% delete the final zeroes
y_hat_ola = y_hat_ola(1:end-wlen);

%% playback: OLA filtered signal y_hat_ola
soundsc(y_hat_ola, Fx);
pause(2);

%% Plot the estimated filter taps against the true RIR
% load the RIR (g.wav)
[g, Fg] = audioread('g.wav');

% produce the plot
time_axis = 0:1/Fg:1/Fg*(M-1);
figure();
subplot(311);
plot(time_axis, g);
title('Real RIR');
xlabel('Time [seconds]');
ylabel('g(n)');
subplot(312);
plot(time_axis, w_o);
title('Wiener-Hopf Taps');
xlabel('Time [seconds]');
ylabel('w_o(n)');
subplot(313);
plot(time_axis, w);
title('Steepest Descent Taps');
xlabel('Time [seconds]');
ylabel('w(n)');

% save current figure as png
saveas(gcf, 'wiener_hopf_estimation_plot.png');

% ----EOF----

