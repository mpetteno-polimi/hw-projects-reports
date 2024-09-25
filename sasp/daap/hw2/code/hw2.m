clc; clear; close all;

%% Your name(s), student ID number(s)
%-------------------------------------------------------------------------%
% Matteo Pettenò, 10868930
% Marco Furio Colombo, 10537094
%-------------------------------------------------------------------------%

% import the array data
load("./array_data_64_mics.mat");

%% Parameters

% sampling frequency
Fs = 8000;

% speed of sound [m/s]
c = 340;

% number of sources
N_src =  2;

% number of microphones
M = 64;

% microphone signal lengthlambda
N = size(y, 2);

% determine the distance between two mics (see anti-aliasing condition)
d = c / Fs;

%% Frequency estimation 
mic_signal = y(1, :);

% plot the spectrum of a microphone signal (use fftshift)
freq_axis = [0:(Fs/N):Fs-(Fs/N)] - Fs/2;
magnitude_spectrum = abs(mic_signal);

figure(1);
plot(freq_axis, fftshift(magnitude_spectrum));
title('Mic signal');
xlabel('Frequency [Hz]');
ylabel('Magnitude');

% estimate of the angular frequency of source signal 1
%[freq_peaks, peaks_index] = maxk(magnitude_spectrum(1:N/2), N_src);
[freq_peaks, peaks_index] = findpeaks(magnitude_spectrum(1:N/2), 'SortStr', 'descend', 'NPeaks', N_src);
f_c1 = (peaks_index(1) - 1)*Fs/N;
w_c1 = 2*pi*f_c1;

% estimate of the angular frequency of source signal 2
f_c2 = (peaks_index(2) - 1)*Fs/N;
w_c2 = 2*pi*f_c2;

% row vector containing the estimated angular frequencies
w_c = [w_c1,w_c2];

%% Delay-and-sum beamformer
% angluar step between candidates for our DOA estimate
candidate_angle_step = deg2rad(1);
% the angles to be examined
candidate_angles = deg2rad(-90):candidate_angle_step:deg2rad(90); 
% the number of the angles to be examined
n_angles = size(candidate_angles, 2);

% pre-allocate the array for the DAS pseudo-spectrum
DAS_pseudo = zeros(n_angles,1);

% pre-allocate the array for the DOA estimates of the two sources 
DAS_DOA = zeros(N_src,1);

% pre-allocate the steering vector for the two sources
a = zeros(M, n_angles, N_src);

% pre-allocate the spatial filter coefficients for the two sources
h_DOA = zeros(M, N_src);

% sample estimate of the covariance matrix of the array data
R = cov(y');

for i = 1:N_src
    
    for j = 1:n_angles
        
        % compute the spatial frequency
        w_s = w_c(i) * d*sin(candidate_angles(j))/c;
        % compute the steering vector for the given source and angle
        a(:,j,i) = exp(-1j*w_s).^(0:M-1).';
        % save the value of the pseudo-spectrum for the given frequency and angle
        DAS_pseudo(j) = a(:,j,i)'*R*a(:,j,i)/M^2;
        
    end

    % remove the spurious imaginary part from the DAS pseudo-spectrum
    DAS_pseudo = real(DAS_pseudo);
    
    % estimate the DOA as the angle for which the DAS pseudo-spectrum has
    % the most prominent peak
    %[DAS_pseudo_max, DAS_pseudo_max_index] = max(DAS_pseudo);
    [DAS_pseudo_max, DAS_pseudo_max_index] = findpeaks(DAS_pseudo, 'SortStr', 'descend', 'NPeaks', 1);
    DAS_DOA(i) = candidate_angles(DAS_pseudo_max_index);
    
    % compute the spatial frequency associated to the DOA
    w_s_DOA = w_c(i) * d*sin(DAS_DOA(i))/c;
    % compute the steering vector associated to the DOA
    a_DOA = exp(-1j*w_s_DOA).^(0:M-1).';
    % compute and save the spatial filter associated to the angle of the DOA
    % for the given source
    h_DOA(:,i) = a_DOA/M; 

    % compute the spatial response using the pre-computed steering
    % vectors and the spatial filter associated to the estimated DOA
    DAS_spatial_response = h_DOA(:,i)' * a(:,:,i);
        
    % plot the DAS pseudo-spectrum for the given source frequency
    figure(i+1)
    subplot(211)
    plot(rad2deg(candidate_angles), DAS_pseudo)
    xlabel('Angle [deg]')
    ylabel('Pseudo-spectrum')
    title("Delay-and-sum beamformer: pseudo-spectrum at " + num2str(w_c(i)/(2*pi)) + " Hz")
    
    % plot the DAS spatial response (beam pattern) for the given source frequency
    figure(i+1)
    subplot(212)
    polarplot(candidate_angles, real(DAS_spatial_response));
    title("Delay-and-sum beamformer: beam pattern at " + num2str(w_c(i)/(2*pi)) + " Hz")
    thetalim([-90, 90])

end

% Display the DOA estimate (in degrees)
DAS_DOA_1_degrees = rad2deg(DAS_DOA(1));
DAS_DOA_2_degrees = rad2deg(DAS_DOA(2));

disp('DELAY-AND-SUM:')
disp(DAS_DOA_1_degrees)
disp(DAS_DOA_2_degrees) 

%% Apply spatial filtering to the array data
mic_time_domain_signal = ifft(mic_signal, N);
s_hat_1 = ifft(h_DOA(:, 1)'*y, N);
s_hat_2 = ifft(h_DOA(:, 2)'*y, N);

%% Play a microphone signal and the two beamformer outputs

% play the mixture (mic signal)
soundsc(mic_time_domain_signal, Fs)
pause
% play the filtered singal (beamformer steered towards source 1)
soundsc(real(s_hat_1), Fs)
pause
% play the filtered singal (beamformer steered towards source 2)
soundsc(real(s_hat_2), Fs)
pause

%% Parametric methods
% eigenvalue decomposition of the covariance matrix R
[Q, R_eigenvalues] = eig(R);

% sort the eigenvalues in descending order
[R_eigenvalues, sorting_index] = sort(diag(R_eigenvalues), 'descend');
    
% permute the columns of the eigenvector matrix accordingly
Q = Q(:, sorting_index);

%% MUSIC
% retain the matrix whose columns span the noise subspace (matrix V)
V = Q(:, N_src+1:M);

% pre-allocate the array for the MUSIC pseudo-spectrum
MUSIC_pseudo = zeros(n_angles,1);
% pre-allocate the array for the DOA estimates of the two sources
MUSIC_DOA = zeros(N_src,1);

for i=1: N_src

    for j=1 : n_angles
        % compute the MUSIC pseudo-spectrum using the precomputed steering vector
        MUSIC_pseudo(j) =  1/(a(:,j,i)'*(V*V')*a(:,j,i));
    end
    
    % remove the spurious imaginary part
    MUSIC_pseudo = real(MUSIC_pseudo);
        
    % estimate the DOA as the angle for which the MUSIC pseudo-spectrum 
    % has the most prominent peak
    %[MUSIC_DOA_pseudo_peak, MUSIC_DOA_pseudo_peak_index] = max(MUSIC_pseudo);
    [MUSIC_DOA_pseudo_peak, MUSIC_DOA_pseudo_peak_index] = findpeaks(MUSIC_pseudo, 'SortStr', 'descend', 'NPeaks', 1);
    MUSIC_DOA(i) = candidate_angles(MUSIC_DOA_pseudo_peak_index);

    % plot the MUSIC pseudo-spectrum for the given source
    figure(4);
    subplot(2,1,i);
    plot(rad2deg(candidate_angles), MUSIC_pseudo);
    xlabel('Angle [deg]');
    ylabel('Pseudo-spectrum');
    title("MUSIC pseudo-spectrum at " + num2str(w_c(i)/(2*pi)) + " Hz");
    
end

% display the MUSIC DOA estimate (in degrees)
MUSIC_DOA_1_degrees = rad2deg(MUSIC_DOA(1));
MUSIC_DOA_2_degrees = rad2deg(MUSIC_DOA(2));

disp('MUSIC:');
disp(MUSIC_DOA_1_degrees);
disp(MUSIC_DOA_2_degrees);

% EOF