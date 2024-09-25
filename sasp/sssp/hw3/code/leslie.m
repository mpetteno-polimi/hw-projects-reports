function [y, y_lpf, y_hpf, y_hp_sdf] = leslie(x, Fs, freq)
% Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);

% allocate output signal
y = zeros(N, 1);

% global modulator parameters
alpha = 0.9;
% tremble spectral delay filter parameter 
Ms_t = 0.2;
Mb_t = -0.75;
N_sdf_t = 4;
% bass spectral delay filter parameter 
Ms_b = 0.04;
Mb_b = -0.92;
N_sdf_b = 3;

% cross-over network design
fc = 800; % cutoff frequency

% TODO: compute the coefficients for the two 4th order butterworth filters
% with cutoff frequency fc
filter_order = 4;
fc_norm = 2 * fc / Fs;
[b_lp, a_lp] = butter(filter_order, fc_norm, "low"); % LPF design
[b_hp, a_hp] = butter(filter_order, fc_norm, "high");  % HPF design
% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state = zeros(length(a_hp)-1, 1);
hpf.in = zeros(length(b_hp), 1);
% lp filter buffers
lpf.state = zeros(length(a_lp)-1, 1);
lpf.in = zeros(length(b_lp), 1); % this could be removed and use only one buffer
% treble sdf filter buffers
sdf_h.state = zeros(N_sdf_t, 1);
sdf_h.in = zeros(N_sdf_t, 1);
% bass sdf filter buffers
sdf_b.state = zeros(N_sdf_b, 1);
sdf_b.in = zeros(N_sdf_b, 1);

% modulators
t = (1:N)/Fs;
m_n_b = sin(2*pi*freq*t);
m_b = Ms_b*m_n_b + Mb_b; % bass modulator
m_n_t = sin(2*pi*(freq + 0.1)*t);
m_t = Ms_t*m_n_t + Mb_t; % treble modulator

% sample processing
for n = 1:N

    % compute crossover network filters output
    lpf.in = circshift(lpf.in, 1);
    lpf.in(1) = x(n);
    y_lpf = (b_lp*lpf.in - a_lp(2:end)*lpf.state) / a_lp(1);
    lpf.state = circshift(lpf.state, 1);
    lpf.state(1) = y_lpf;

    hpf.in = circshift(hpf.in , 1);
    hpf.in(1) = x(n);
    y_hpf = (b_hp*hpf.in - a_hp(2:end)*hpf.state) / a_hp(1);
    hpf.state = circshift(hpf.state, 1);
    hpf.state(1) = y_hpf;

    % compute bass SDF output
    y_lp_sdf = sdf_b.in(end);
    sdf_b.in = circshift(sdf_b.in, 1);
    sdf_b.in(1) = y_lpf;
    for i = 1:N_sdf_b
        y_lp_sdf = y_lp_sdf + nchoosek(N_sdf_b, i)*m_b(n)^i*(sdf_b.in(end-i+1)-sdf_b.state(i));
    end
    sdf_b.state = circshift(sdf_b.state, 1);
    sdf_b.state(1) = y_lp_sdf;

    % compute treble SDF output
    y_hp_sdf = sdf_h.in(end);
    sdf_h.in = circshift(sdf_h.in, 1);
    sdf_h.in(1) = y_hpf;
    for i = 1:N_sdf_t 
        y_hp_sdf = y_hp_sdf + nchoosek(N_sdf_t, i)*m_t(n)^i*(sdf_h.in(end-i+1)-sdf_h.state(i));
    end
    sdf_h.state = circshift(sdf_h.state, 1);
    sdf_h.state(1) = y_hp_sdf;

    % implement AM modulation blocks
    y_lp_am = (1 + alpha*m_b(n))*y_lp_sdf;
    y_hp_am = (1 + alpha*m_t(n))*y_hp_sdf;

    y(n) = y_lp_am + y_hp_am;

end

end

