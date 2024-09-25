%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2022                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
Fs = 48000;                         % sampling frequency
Ts = 1/Fs;                          % temporal resolution (time step)
duration = 8;                       % total length is 8s
N = duration/Ts;                    % number of temporal samples
t_axis = linspace(0, duration, N);  % time axis
f_axis = 0:Fs/N:Fs-Fs/N;            % frequency axis

% Fundamental note
f_1 = 65.4;

% Boundary
nu_l = 1e20;                        % left end normalized impedance
nu_b = 1e3;                         % bridge normalized impedance

% String parameters
M_s = 35e-3;                        % string mass
L = 1.92;                           % string length
rho = M_s/L;                        % string linear density 
T_e = rho*(2*L*f_1)^2;              % string tension
b_1 = 0.5;                          % air damping coefficient
b_2 = 6.25e-9;                      % string internal friction coefficient
eps = 7.5e-6;                       % string stiffness parameter

% Spatial sampling parameters
% Aliasing condition
% Number of maximum spatial steps
gamma = Fs/(2*f_1);
M_max = sqrt((-1+sqrt(1+16*eps*gamma^2))/(8*eps));
% Integer values
M = floor(M_max);                   % number of string segment
% Spatial sampling
Xs = L/M;                           % spatial resolution (spatial step)

% FD parameters
c = sqrt(T_e/rho);                  % string transverse wave velocity
lambda = c*Ts/Xs;                   % Courant number < 1
mu = eps^2/(c^2*Xs^2);
v = (2*b_2*Ts)/(Xs^2);

% Hammer parametersn
M_h = 4.9e-3;                       % hammer mass
b_h = 1e-4;                         % fluid damping coefficient
p = 2.3;                            % stiffness exponent
K = 4e8;                            % stiffness coefficient
w = 0.2;                            % width of the hammer spatial window g
a_h = 0.12;                         % relative striking position (L/x_0)

% Hammer contact window definition
x_0 = a_h*L;                        % excitation point
m_0 = round(x_0/Xs);                % discretize x_0 to closest integer
L_g = 2*floor(w/(2*Xs))+1;          % discretize w to closest odd integer
g = zeros(1, M);
g((m_0-(L_g-1)/2):m_0+(L_g-1)/2) = hann(L_g);

% PDE Coefficients
a_1 = (-lambda^2*mu)/(1+b_1*Ts);
a_2 = (lambda^2+4*lambda^2*mu+v)/(1+b_1*Ts);
a_3 = (2-2*lambda^2-6*lambda^2*mu-2*v)/(1+b_1*Ts);
a_4 = (-1+b_1*Ts+2*v)/(1+b_1*Ts);
a_5 = -v/(1+b_1*Ts);
a_f = (Ts^2/rho)/(1+b_1*Ts);

% Bridge boundary coefficients
b_r1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*Ts+nu_b*lambda);
b_r2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*Ts+nu_b*lambda);
b_r3 = (-2*lambda^2*mu)/(1+b_1*Ts+nu_b*lambda);
b_r4 = (-1+b_1*Ts+nu_b*lambda)/(1+b_1*Ts+nu_b*lambda);
b_rf = (Ts^2/rho)/(1+b_1*Ts+nu_b*lambda);

% Left hand (hinged string end) boundary coefficients
b_l1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*Ts+nu_l*lambda);
b_l2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*Ts+nu_l*lambda);
b_l3 = (-2*lambda^2*mu)/(1+b_1*Ts+nu_l*lambda);
b_l4 = (-1+b_1*Ts+nu_l*lambda)/(1+b_1*Ts+nu_l*lambda);
b_lf = (Ts^2/rho)/(1+b_1*Ts+nu_l*lambda);

% Hammer felt parameters
d_1 = 2/(1+b_h*Ts/(2*M_h));
d_2 = (-1+b_h*Ts/(2*M_h))/(1+b_h*Ts/(2*M_h));
d_f = (-Ts^2/M_h)/(1+b_h*Ts/(2*M_h));

%% Computation of the FD scheme

% Initialization
y = zeros(M, N);                        % initial string displacement
soundwave = zeros(1, N);                % initial string velocity
V_ho = 2.5;                             % initial hammer velocity
eta = zeros(1, N); eta(2) = V_ho*Ts;    % initial hammer displacement
F_h = zeros(M, N);                      % initial hammer force
X_av = 12;
hammer_on = true;

% Computation loop
for n = 2:N-1
    % The interaction process ends when the displacement of the hammer head 
    % becomes less than the displacement of the string at the center of the
    % contact segment (x0)
    if (hammer_on && eta(n) < y(m_0,n))
        hammer_on = false;
        F_h_n = 0;
    end
    
    if (hammer_on)
        F_h_n = K*abs(eta(n)-y(m_0,n))^p;
        % Compute hammer displacement for next step
        eta(n+1) = d_1*eta(n)+d_2*eta(n-1)+d_f*F_h_n;
    end
    
    % Compute spatial samples for current time
    for m = 1:M
        F_h(m,n) = F_h_n*g(m);
        if (m == 1)
            y(m,n+1) = b_l1*y(m,n) + b_l2*y(m+1,n) + b_l3*y(m+2,n) + ...
                b_l4*y(m,n-1) + b_lf*F_h(m,n);
        elseif (m == M)
            y(m,n+1) = b_r1*y(m,n) + b_r2*y(m-1,n) + b_r3*y(m-2,n) + ...
                b_r4*y(m,n-1) + b_rf*F_h(m,n);
        elseif (m == 2)
            y(m,n+1) = a_1*(y(m+2,n)-y(m,n)+2*y(m-1,n)) + ...
                a_2*(y(m+1,n)+y(m-1,n)) + a_3*y(m,n) + a_4*y(m,n-1) + ...
                a_5*(y(m+1,n-1)+y(m-1,n-1)) + a_f*F_h(m,n);
        elseif (m == M-1)
            y(m,n+1) = a_1*(2*y(m+1,n)-y(m,n)+y(m-2,n)) + ...
                a_2*(y(m+1,n)+y(m-1,n)) + a_3*y(m,n) + a_4*y(m,n-1) + ...
                a_5*(y(m+1,n-1)+y(m-1,n-1)) + a_f*F_h(m,n);
        else
            y(m,n+1) = a_1*(y(m+2,n)+y(m-2,n)) + a_2*(y(m+1,n)+y(m-1,n)) + ...
                a_3*y(m,n) + a_4*y(m,n-1) + a_5*(y(m+1,n-1)+y(m-1,n-1)) + ...
                a_f*F_h(m,n);
        end
    end
    
    % Average the sound over X_av spatial samples centered in x_0
    soundwave(n) = sum(y((m_0-X_av/2):m_0+X_av/2, n))/X_av;
end

%% Plot the displacement in time
figure();
plot(t_axis, y(m_0, :));
ylabel('Displacement in x_0 [m]');
xlabel('Time [s]');

%% Plot the synthesized signal play it and save it on the disk
figure();
plot(t_axis, soundwave);
ylabel('Synthesized signal');
xlabel('Time [s]');

figure();
plot(f_axis, db(abs(fft(soundwave))));
ylabel('Magnitude [dB]');
xlabel('Frequency [Hz]');
xlim([0, 1500])

% Play the sound
synth_piano = soundwave./max(abs(soundwave));
sound(synth_piano, Fs);

% Save on disk
disp('Save file on disk...');
id_number = '10868930';
surname = 'petteno';
fileName = sprintf('%s_%s_piano.wav', id_number, surname);
audiowrite(fileName, synth_piano, Fs);
